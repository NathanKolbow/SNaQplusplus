#!/usr/bin/env python3
"""CSNaQ benchmark framework.

Runs fixed benchmark suites on a `snaq` binary (or on SNaQ.jl, the baseline),
appends the results to a running history (benchmarks/history/results.jsonl),
and draws comparison figures and a report (benchmarks/reports/).

Only the Python standard library is needed. Running the SNaQ.jl baseline
additionally needs Julia with the environment of tools/julia
(see tools/julia/setup.jl). See benchmarks/README.md for usage.
"""

import argparse
import datetime as _dt
import json
import math
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = os.path.join(ROOT, "benchmarks")
DATASETS = os.path.join(BENCH, "datasets")
HISTORY_DIR = os.path.join(BENCH, "history")
HISTORY = os.path.join(HISTORY_DIR, "results.jsonl")
OUTPUTS = os.path.join(HISTORY_DIR, "outputs")
REPORTS = os.path.join(BENCH, "reports")
JULIA_DIR = os.path.join(BENCH, "julia")
JULIA_PROJECT = os.path.join(ROOT, "tools", "julia")

DEFAULT_SNAQ = os.path.join(ROOT, "build", "snaq")
DEFAULT_EVAL = os.path.join(ROOT, "build", "csnaq_eval")
DEFAULT_FIT = os.path.join(ROOT, "build", "csnaq_bench_fit")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def load_suites():
    with open(os.path.join(BENCH, "suites.json")) as f:
        return json.load(f)


def exe(path):
    """Path of an executable, adding .exe (and Release/) on Windows if needed."""
    if os.path.isfile(path):
        return path
    for cand in (path + ".exe", os.path.join(os.path.dirname(path), "Release", os.path.basename(path) + ".exe")):
        if os.path.isfile(cand):
            return cand
    return path


def git_info():
    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
        except Exception:
            return ""
    commit = git("rev-parse", "--short", "HEAD")
    status = git("status", "--porcelain", "--", ".", ":(exclude)benchmarks/history", ":(exclude)benchmarks/reports")
    return commit, bool(status)


def host_info():
    cpu = platform.processor()
    try:
        if sys.platform == "darwin":
            cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True,
                                 text=True).stdout.strip() or cpu
        elif sys.platform.startswith("linux"):
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        cpu = line.split(":", 1)[1].strip()
                        break
    except Exception:
        pass
    return {"os": f"{platform.system()} {platform.release()}", "machine": platform.machine(), "cpu": cpu,
            "cpu_count": os.cpu_count(), "hostname": platform.node()}


def run_measured(cmd, cwd=None):
    """Run a command; return (returncode, stdout, stderr, wall seconds, peak RSS in MB or None)."""
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = [], []
    readers = [threading.Thread(target=lambda: out.append(p.stdout.read())),
               threading.Thread(target=lambda: err.append(p.stderr.read()))]
    for r in readers:
        r.start()
    peak = None
    if hasattr(os, "wait4"):
        _, status, ru = os.wait4(p.pid, 0)
        rc = os.waitstatus_to_exitcode(status)
        p.returncode = rc
        # ru_maxrss: bytes on macOS, kilobytes on Linux
        peak = ru.ru_maxrss / 2**20 if sys.platform == "darwin" else ru.ru_maxrss / 2**10
    else:
        rc = p.wait()
    wall = time.perf_counter() - t0
    for r in readers:
        r.join()
    return rc, "".join(out), "".join(err), wall, peak


def parse_out_file(path):
    """Best network and scores of all runs from a SNaQ .out file."""
    with open(path) as f:
        lines = f.read().splitlines()
    best = lines[0].split(" SNaQ score = ")[0]
    best_score = float(lines[0].split(" SNaQ score = ")[1])
    scores = [float(m.group(1)) for l in lines if (m := re.search(r", with SNaQ score (\S+)$", l))]
    return best, best_score, scores


def evaluate(eval_bin, estimated, truth, genetrees, threads):
    rc, out, err, _, _ = run_measured([exe(eval_bin), estimated, truth, genetrees, str(threads)])
    if rc != 0:
        print(f"    warning: accuracy evaluation failed: {err.strip()}")
        return {}
    return json.loads(out)


def dataset_files(name):
    d = os.path.join(DATASETS, name)
    return {k: os.path.join(d, f) for k, f in
            (("truth", "truth.net"), ("genetrees", "genetrees.tre"), ("start", "start.tre"))}


def dataset_ntaxa(files):
    """Number of taxa of a dataset (leaf names of its starting tree)."""
    with open(files["start"]) as f:
        newick = f.readline()
    return len(set(re.findall(r"[(,]\s*([^(),:;#\s\[]+)", newick)))


def slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", s).strip("-") or "run"


def append_record(rec):
    os.makedirs(HISTORY_DIR, exist_ok=True)
    with open(HISTORY, "a") as f:
        f.write(json.dumps(rec, sort_keys=True) + "\n")


def load_history():
    if not os.path.exists(HISTORY):
        return []
    with open(HISTORY) as f:
        return [json.loads(l) for l in f if l.strip()]


def select_cases(suite, names):
    cases = suite["cases"]
    if names:
        wanted = set(names.split(","))
        cases = [c for c in cases if c["name"] in wanted]
        missing = wanted - {c["name"] for c in cases}
        if missing:
            sys.exit(f"unknown case(s): {', '.join(sorted(missing))}")
    return cases


# --------------------------------------------------------------------------
# running benchmarks
# --------------------------------------------------------------------------

def run_suite(args, implementation):
    suites = load_suites()
    if args.suite not in suites:
        sys.exit(f"unknown suite {args.suite}; available: {', '.join(suites)}")
    suite = suites[args.suite]
    cases = select_cases(suite, args.cases)
    commit, dirty = git_info()
    label = args.label or (("SNaQ.jl" if implementation == "julia" else "cpp") + "@" + commit + ("+dirty" if dirty else ""))
    now = _dt.datetime.now().astimezone()
    run_id = now.strftime("%Y%m%dT%H%M%S") + "_" + slug(label)
    outdir = os.path.join(OUTPUTS, run_id)
    os.makedirs(outdir, exist_ok=True)
    common = {"run_id": run_id, "label": label, "implementation": implementation, "suite": args.suite,
              "kind": suite["kind"], "timestamp": now.isoformat(timespec="seconds"), "git_commit": commit,
              "git_dirty": dirty, "threads": args.threads, "notes": args.notes or "", "host": host_info()}
    if implementation == "cpp":
        binary = exe(args.binary if suite["kind"] == "snaq" else args.fit_binary)
        if not os.path.isfile(binary):
            sys.exit(f"binary not found: {binary} (build first: see benchmarks/README.md)")
        common["binary"] = os.path.relpath(binary, ROOT) if binary.startswith(ROOT) else binary
        common["extra_args"] = args.extra_args or ""
    print(f"run {run_id}: suite '{args.suite}', {len(cases)} case(s), label '{label}'")
    for case in cases:
        files = dataset_files(case["dataset"])
        case = dict(case)
        if suite["kind"] == "snaq" and "nquartets" in suite and "propQuartets" not in case:
            # about `nquartets` quartets per case: min(1, nquartets / binomial(ntaxa, 4))
            ntaxa = dataset_ntaxa(files)
            # (the -1e-9 keeps ceil(binomial(ntaxa,4) * prop) at exactly `nquartets`
            # despite floating-point rounding; both implementations sample that many)
            prop = min(1.0, (suite["nquartets"] - 1e-9) / math.comb(ntaxa, 4))
            case.setdefault("propQuartets", prop)
            case.setdefault("propQuartetsFinal", prop)
        print(f"  {case['name']} ...", flush=True)
        rec = dict(common)
        rec["case"] = case["name"]
        rec["params"] = case
        try:
            if suite["kind"] == "snaq":
                rec["metrics"] = run_snaq_case(args, implementation, case, files, outdir)
            else:
                rec["metrics"] = run_fit_case(args, implementation, case, files)
        except RuntimeError as e:
            print(f"    FAILED: {e}")
            rec["metrics"] = {"failed": True, "error": str(e)[-2000:]}
        append_record(rec)
        m = rec["metrics"]
        if not m.get("failed"):
            if suite["kind"] == "snaq":
                print(f"    time {m['time_s']:.2f} s, peak memory {fmt_mb(m.get('peak_rss_mb'))}, "
                      f"best score {m.get('best_score')}, HWCD to truth {m.get('hwcd')}")
            else:
                print(f"    fit {m['fit_s']:.4f} s, equations {m['equations_s']:.4f} s, score {m['score']}, "
                      f"peak memory {fmt_mb(m.get('peak_rss_mb'))}")
    print(f"results appended to {os.path.relpath(HISTORY, ROOT)}; outputs in {os.path.relpath(outdir, ROOT)}")


def fmt_mb(x):
    return "n/a" if x is None else f"{x:.1f} MB"


def run_snaq_case(args, implementation, case, files, outdir):
    pq = case.get("propQuartets", 1.0)
    pqf = case.get("propQuartetsFinal", 1.0)
    lazy = bool(case.get("lazy", False))
    seed = case.get("seed", 42)
    tmp = tempfile.mkdtemp(prefix="csnaq_bench_")
    outroot = os.path.join(tmp, "snaq")
    try:
        if implementation == "cpp":
            cmd = [exe(args.binary), "-i", files["genetrees"], "--start-tree", files["start"], "-o", outroot,
                   "--hmax", str(case["hmax"]), "--runs", str(case["runs"]), "--Nfail", str(case["Nfail"]),
                   "--seed", str(seed), "--threads", str(args.threads), "--propQuartets", repr(pq),
                   "--propQuartetsFinal", repr(pqf)]
            if lazy:
                cmd.append("--lazy")
            cmd += shlex.split(args.extra_args or "")
            rc, out, err, wall, peak = run_measured(cmd)
            if rc != 0:
                raise RuntimeError(err.strip() or out.strip())
            metrics = {"time_s": wall, "process_wall_s": wall, "peak_rss_mb": peak}
        else:
            cmd = [args.julia, f"-t{args.threads}", f"--project={JULIA_PROJECT}",
                   os.path.join(JULIA_DIR, "run_snaq_case.jl"), files["genetrees"], files["start"],
                   str(case["hmax"]), str(case["runs"]), str(case["Nfail"]), str(seed), repr(pq), repr(pqf),
                   "true" if lazy else "false", outroot]
            rc, out, err, wall, peak = run_measured(cmd)
            m = re.search(r"RESULT (\{.*\})", out)
            if rc != 0 or not m:
                raise RuntimeError(err.strip()[-2000:] or out.strip()[-2000:])
            r = json.loads(m.group(1))
            metrics = {"time_s": r["snaq_s"] + r["obscf_s"], "snaq_s": r["snaq_s"], "obscf_s": r["obscf_s"],
                       "process_wall_s": wall, "peak_rss_mb": peak if peak is not None else r["peak_rss_mb"],
                       "base_rss_mb": r["base_rss_mb"]}
        best, best_score, scores = parse_out_file(outroot + ".out")
        metrics.update({"best_score": best_score, "run_scores": scores,
                        "median_run_score": statistics.median(scores) if scores else None,
                        "best_network": best})
        if args.eval_binary and os.path.isfile(exe(args.eval_binary)):
            ev = evaluate(args.eval_binary, outroot + ".out", files["truth"], files["genetrees"], args.threads)
            metrics.update({"hwcd": ev.get("hwcd"), "score_all_quartets": ev.get("score_estimated"),
                            "score_truth_all_quartets": ev.get("score_truth"), "nhybrids": ev.get("nhybrids")})
        dest = os.path.join(outdir, case["name"])
        for ext in (".out", ".networks"):
            if os.path.exists(outroot + ext):
                shutil.copy(outroot + ext, dest + ext)
        return metrics
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_fit_case(args, implementation, case, files):
    repeats = str(args.repeats)
    if implementation == "cpp":
        # extra arguments go after the repeats (e.g. --extra-args newton: the optimization method)
        cmd = [exe(args.fit_binary), files["truth"], files["genetrees"], str(case["maxeval"]), str(args.threads),
               repeats] + shlex.split(args.extra_args or "")
        rc, out, err, wall, peak = run_measured(cmd)
        if rc != 0:
            raise RuntimeError(err.strip() or out.strip())
        v = out.split()
        return {"fit_s": float(v[5]), "equations_s": float(v[6]), "score": float(v[4]), "nparams": int(v[3]),
                "nquartets": int(v[2]), "peak_rss_mb": peak, "process_wall_s": wall}
    cmd = [args.julia, f"-t{args.threads}", f"--project={JULIA_PROJECT}", os.path.join(JULIA_DIR, "run_fit_case.jl"),
           files["truth"], files["genetrees"], str(case["maxeval"]), repeats]
    rc, out, err, wall, peak = run_measured(cmd)
    m = re.search(r"RESULT (\{.*\})", out)
    if rc != 0 or not m:
        raise RuntimeError(err.strip()[-2000:] or out.strip()[-2000:])
    r = json.loads(m.group(1))
    r["process_wall_s"] = wall
    if peak is not None:
        r["peak_rss_mb"] = peak
    return r


# --------------------------------------------------------------------------
# listing
# --------------------------------------------------------------------------

def list_runs(args):
    recs = load_history()
    runs = {}
    for r in recs:
        runs.setdefault(r["run_id"], []).append(r)
    if not runs:
        print("no benchmark results yet")
        return
    print(f"{'run id':40} {'label':28} {'suite':10} {'cases':>5} {'failed':>6}  commit")
    for rid, rs in runs.items():
        r0 = rs[0]
        nfail = sum(1 for r in rs if r["metrics"].get("failed"))
        print(f"{rid:40} {r0['label'][:28]:28} {r0['suite']:10} {len(rs):5d} {nfail:6d}  "
              f"{r0['git_commit']}{'+dirty' if r0['git_dirty'] else ''}")


# --------------------------------------------------------------------------
# SVG figures (no external dependencies)
# --------------------------------------------------------------------------

# Reference palette (light): categorical slots in fixed order, and chart chrome.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
FONT = "font-family=\"-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif\""


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;"))


def fmt_num(x):
    if x is None:
        return "n/a"
    ax = abs(x)
    if ax == 0:
        return "0"
    if ax >= 1000:
        return f"{x:,.0f}"
    if ax >= 100:
        return f"{x:.0f}"
    if ax >= 10:
        return f"{x:.1f}"
    if ax >= 1:
        return f"{x:.2f}"
    if ax >= 0.01:
        return f"{x:.3f}"
    return f"{x:.2g}"


def log_ticks(lo, hi):
    a, b = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
    ticks = []
    for e in range(a, b + 1):
        for m in (1, 2, 5):
            v = m * 10 ** e
            if lo <= v <= hi:
                ticks.append(v)
    if len(ticks) > 9:
        ticks = [t for t in ticks if str(t).lstrip("0.").startswith("1") or math.log10(t).is_integer()]
    return ticks


def lin_ticks(lo, hi, n=5):
    span = hi - lo if hi > lo else 1.0
    step = 10 ** math.floor(math.log10(span / n))
    for m in (1, 2, 5, 10):
        if span / (step * m) <= n:
            step *= m
            break
    t = math.ceil(lo / step) * step
    ticks = []
    while t <= hi + 1e-12:
        ticks.append(round(t, 12))
        t += step
    return ticks


def legend_svg(labels, x, y, maxwidth):
    parts, cx, cy = [], x, y
    for i, lab in enumerate(labels):
        w = 22 + 7 * len(lab)
        if cx + w > x + maxwidth and cx > x:
            cx, cy = x, cy + 20
        parts.append(f'<rect x="{cx}" y="{cy - 9}" width="12" height="12" rx="3" fill="{SERIES[i % 8]}"/>')
        parts.append(f'<text x="{cx + 18}" y="{cy + 1}" font-size="12" fill="{INK2}">{esc(lab)}</text>')
        cx += w + 14
    return "\n".join(parts), cy - y + 20


def grouped_bar_svg(title, subtitle, categories, series_labels, values, unit, log=False, reference=None,
                    reference_label=None):
    """Horizontal grouped bars: values[s][c] for series s and category c (None = missing)."""
    W = 860
    left, right = 170, 90
    bar, gap, groupgap = 12, 2, 14
    plotw = W - left - right
    legend, legend_h = legend_svg(series_labels, left, 70, plotw) if len(series_labels) > 1 else ("", 0)
    top = 70 + legend_h + 10
    nser = len(series_labels)
    group_h = nser * bar + (nser - 1) * gap
    H = top + len(categories) * (group_h + groupgap) + 50
    vals = [v for row in values for v in row if v is not None and (v > 0 or not log)]
    if reference is not None:
        vals.append(reference)
    if not vals:
        return None
    if log:
        lo, hi = min(vals) / 1.5, max(vals) * 1.5
        ticks = log_ticks(lo, hi)
        lo, hi = min([lo] + ticks), max([hi] + ticks)

        def X(v):
            return left + plotw * (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
        base = lo
    else:
        lo, hi = min(0.0, min(vals)), max(0.0, max(vals))
        pad = 0.08 * ((hi - lo) or 1.0)
        lo, hi = (lo - pad if lo < 0 else lo), hi + pad
        ticks = lin_ticks(lo, hi)
        lo, hi = min([lo] + ticks), max([hi] + ticks)

        def X(v):
            return left + plotw * (v - lo) / (hi - lo)
        base = 0.0
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" {FONT}>',
           f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
           f'<text x="20" y="30" font-size="17" font-weight="600" fill="{INK}">{esc(title)}</text>',
           f'<text x="20" y="52" font-size="12.5" fill="{INK2}">{esc(subtitle)}</text>', legend]
    axis_y = H - 40
    for t in ticks:
        x = X(t)
        out.append(f'<line x1="{x:.1f}" y1="{top - 6}" x2="{x:.1f}" y2="{axis_y}" stroke="{GRID}" stroke-width="1"/>')
        out.append(f'<text x="{x:.1f}" y="{axis_y + 16}" font-size="11" fill="{MUTED}" text-anchor="middle">'
                   f'{esc(fmt_num(t))}</text>')
    out.append(f'<text x="{left + plotw / 2:.1f}" y="{axis_y + 34}" font-size="11.5" fill="{INK2}" '
               f'text-anchor="middle">{esc(unit)}{" (log scale)" if log else ""}</text>')
    x0 = X(base)
    out.append(f'<line x1="{x0:.1f}" y1="{top - 6}" x2="{x0:.1f}" y2="{axis_y}" stroke="{AXIS}" stroke-width="1"/>')
    if reference is not None:
        xr = X(reference)
        out.append(f'<line x1="{xr:.1f}" y1="{top - 6}" x2="{xr:.1f}" y2="{axis_y}" stroke="{INK2}" '
                   f'stroke-width="1" stroke-dasharray="4 3"/>')
        if reference_label:
            out.append(f'<text x="{xr + 4:.1f}" y="{top - 10}" font-size="11" fill="{INK2}">{esc(reference_label)}</text>')
    for ci, cat in enumerate(categories):
        gy = top + ci * (group_h + groupgap)
        out.append(f'<text x="{left - 10}" y="{gy + group_h / 2 + 4:.1f}" font-size="12" fill="{INK}" '
                   f'text-anchor="end">{esc(cat)}</text>')
        for si in range(nser):
            v = values[si][ci]
            y = gy + si * (bar + gap)
            if v is None or (log and v <= 0):
                out.append(f'<text x="{x0 + 4:.1f}" y="{y + bar - 2}" font-size="10.5" fill="{MUTED}">n/a</text>')
                continue
            xv = X(v)
            xa, xb = (x0, xv) if xv >= x0 else (xv, x0)
            w = max(xb - xa, 1.0)
            r = min(4, w / 2)
            out.append(f'<g><title>{esc(series_labels[si])} · {esc(cat)}: {esc(fmt_num(v))} {esc(unit)}</title>'
                       f'<rect x="{xa:.1f}" y="{y}" width="{w:.1f}" height="{bar}" rx="{r:.1f}" '
                       f'fill="{SERIES[si % 8]}"/></g>')
            tx, anchor = (xb + 5, "start") if xv >= x0 else (xa - 5, "end")
            out.append(f'<text x="{tx:.1f}" y="{y + bar - 2}" font-size="10.5" fill="{INK2}" '
                       f'text-anchor="{anchor}">{esc(fmt_num(v))}</text>')
    out.append("</svg>")
    return "\n".join(out)


def line_svg(title, subtitle, xlabels, series_labels, values, unit, log=True, reference=None, reference_label=None):
    """Lines over a categorical x axis: values[s][i] (None = missing)."""
    W, left, right = 860, 70, 30
    plotw = W - left - right
    legend, legend_h = legend_svg(series_labels, left, 70, plotw) if len(series_labels) > 1 else ("", 0)
    top = 70 + legend_h + 14
    ploth = 260
    H = top + ploth + 90
    vals = [v for row in values for v in row if v is not None and (v > 0 or not log)]
    if reference is not None:
        vals.append(reference)
    if not vals:
        return None
    if log:
        lo, hi = min(vals) / 1.5, max(vals) * 1.5
        ticks = log_ticks(lo, hi)
        lo, hi = min([lo] + ticks), max([hi] + ticks)

        def Y(v):
            return top + ploth * (1 - (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)))
    else:
        lo, hi = min(vals), max(vals)
        ticks = lin_ticks(lo, hi)
        lo, hi = min([lo] + ticks), max([hi] + ticks)

        def Y(v):
            return top + ploth * (1 - (v - lo) / ((hi - lo) or 1))
    n = len(xlabels)

    def X(i):
        return left + (plotw * (i + 0.5) / n)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" {FONT}>',
           f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
           f'<text x="20" y="30" font-size="17" font-weight="600" fill="{INK}">{esc(title)}</text>',
           f'<text x="20" y="52" font-size="12.5" fill="{INK2}">{esc(subtitle)}</text>', legend]
    for t in ticks:
        y = Y(t)
        out.append(f'<line x1="{left}" y1="{y:.1f}" x2="{W - right}" y2="{y:.1f}" stroke="{GRID}" stroke-width="1"/>')
        out.append(f'<text x="{left - 8}" y="{y + 4:.1f}" font-size="11" fill="{MUTED}" text-anchor="end">'
                   f'{esc(fmt_num(t))}</text>')
    out.append(f'<text x="16" y="{top + ploth / 2:.1f}" font-size="11.5" fill="{INK2}" text-anchor="middle" '
               f'transform="rotate(-90 16 {top + ploth / 2:.1f})">{esc(unit)}{" (log)" if log else ""}</text>')
    if reference is not None:
        yr = Y(reference)
        out.append(f'<line x1="{left}" y1="{yr:.1f}" x2="{W - right}" y2="{yr:.1f}" stroke="{INK2}" '
                   f'stroke-width="1" stroke-dasharray="4 3"/>')
        if reference_label:
            out.append(f'<text x="{left + 6}" y="{yr - 5:.1f}" font-size="11" fill="{INK2}">'
                       f'{esc(reference_label)}</text>')
    for i, lab in enumerate(xlabels):
        x = X(i)
        out.append(f'<text x="{x:.1f}" y="{top + ploth + 16}" font-size="10.5" fill="{MUTED}" text-anchor="end" '
                   f'transform="rotate(-30 {x:.1f} {top + ploth + 16})">{esc(lab[:32])}</text>')
    for si, row in enumerate(values):
        pts = [(X(i), Y(v), v, xlabels[i]) for i, v in enumerate(row) if v is not None and (v > 0 or not log)]
        if len(pts) > 1:
            d = " ".join(("M" if k == 0 else "L") + f"{x:.1f},{y:.1f}" for k, (x, y, _, _) in enumerate(pts))
            out.append(f'<path d="{d}" fill="none" stroke="{SERIES[si % 8]}" stroke-width="2" stroke-linejoin="round"/>')
        for x, y, v, xl in pts:
            out.append(f'<g><title>{esc(series_labels[si])} · {esc(xl)}: {esc(fmt_num(v))} {esc(unit)}</title>'
                       f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{SERIES[si % 8]}" stroke="{SURFACE}" '
                       f'stroke-width="2"/></g>')
        if pts:  # direct label at the last point
            x, y, v, _ = pts[-1]
            out.append(f'<text x="{x + 8:.1f}" y="{y + 4:.1f}" font-size="10.5" fill="{INK2}">{esc(fmt_num(v))}</text>')
    out.append("</svg>")
    return "\n".join(out)


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def latest_runs_by_label(recs, suite):
    """label -> {case -> record} using each label's most recent run of the suite."""
    by_label = {}
    for r in recs:
        if r["suite"] != suite:
            continue
        cur = by_label.get(r["label"])
        if cur is None or r["run_id"] > cur:
            by_label[r["label"]] = r["run_id"]
    out = {}
    for label, rid in by_label.items():
        out[label] = {r["case"]: r for r in recs if r["run_id"] == rid}
    return out


def ordered_labels(recs, suite, baseline, chosen, maxseries=4):
    first_seen = {}
    for r in sorted(recs, key=lambda r: r["run_id"]):
        if r["suite"] == suite:
            first_seen.setdefault(r["label"], r["run_id"])
    labels = sorted(first_seen, key=lambda l: first_seen[l])
    if chosen:
        return [l for l in chosen if l in first_seen]
    rest = [l for l in labels if l != baseline]
    keep = rest[-(maxseries - 1):] if baseline in first_seen else rest[-maxseries:]
    return ([baseline] if baseline in first_seen else []) + keep


RESCORED_KEYS = ("hwcd", "score_all_quartets", "score_truth_all_quartets", "nhybrids")


def metric(rec, key):
    """A metric of a record. Accuracy metrics come from `bench.py rescore` when it was run,
    so that every version is evaluated by the same csnaq_eval."""
    if rec is None or rec["metrics"].get("failed"):
        return None
    m = rec["metrics"]
    if key == "score_gap":
        a, b = metric(rec, "score_all_quartets"), metric(rec, "score_truth_all_quartets")
        return None if a is None or b is None else a - b
    if key in RESCORED_KEYS and key in m.get("rescored", {}):
        return m["rescored"][key]
    return m.get(key)


def rescore(args):
    """Recompute the accuracy metrics of every recorded snaq run from its saved .out file."""
    recs = load_history()
    eval_bin = exe(args.eval_binary)
    if not os.path.isfile(eval_bin):
        sys.exit(f"csnaq_eval not found: {eval_bin} (build the tools first)")
    commit, dirty = git_info()
    suites = load_suites()
    n = 0
    done = {}
    for r in recs:
        if r["kind"] != "snaq" or r["metrics"].get("failed"):
            continue
        out = os.path.join(OUTPUTS, r["run_id"], r["case"] + ".out")
        if not os.path.isfile(out):
            print(f"  {r['run_id']} {r['case']}: no saved output, skipped")
            continue
        files = dataset_files(r["params"]["dataset"])
        ev = evaluate(eval_bin, out, files["truth"], files["genetrees"], args.threads)
        if not ev:
            continue
        done[(r["run_id"], r["case"])] = {"hwcd": ev.get("hwcd"), "score_all_quartets": ev.get("score_estimated"),
                                           "score_truth_all_quartets": ev.get("score_truth"),
                                           "nhybrids": ev.get("nhybrids"),
                                           "eval_commit": commit + ("+dirty" if dirty else "")}
        n += 1
    # read the history again: records may have been appended in the meantime
    recs = load_history()
    for r in recs:
        if (r["run_id"], r["case"]) in done:
            r["metrics"]["rescored"] = done[(r["run_id"], r["case"])]
    tmp = HISTORY + ".tmp"
    with open(tmp, "w") as f:
        for r in recs:
            f.write(json.dumps(r, sort_keys=True) + "\n")
    os.replace(tmp, HISTORY)
    print(f"rescored {n} record(s) with {os.path.relpath(eval_bin, ROOT) if eval_bin.startswith(ROOT) else eval_bin}")


def geomean(xs):
    xs = [x for x in xs if x is not None and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def write(path, text):
    if text is None:
        return False
    with open(path, "w") as f:
        f.write(text)
    return True


def report(args):
    recs = load_history()
    if not recs:
        sys.exit("no benchmark results yet: run `bench.py run` first")
    suites = load_suites()
    os.makedirs(REPORTS, exist_ok=True)
    chosen = args.labels.split(",") if args.labels else None
    md = ["# CSNaQ benchmark report", "",
          f"Generated {_dt.datetime.now().astimezone().isoformat(timespec='seconds')} by `benchmarks/bench.py report` "
          f"from `benchmarks/history/results.jsonl`. Baseline: **{args.baseline}**.", "",
          "Time is the wall time of the whole `snaq` process for C++. For SNaQ.jl it is the time to build "
          "the DataCF plus the `snaq!` call, after a warm-up call, so Julia's start-up and compilation are "
          "not included. Peak memory is the peak resident set size of the process; for SNaQ.jl it includes "
          "the Julia runtime and loaded packages (see `base_rss_mb` in the history). HWCD is the hardwired "
          "cluster distance between the best network and the true network (0: same topology). The score on "
          "all quartets is the SNaQ score of the best network recomputed by `csnaq_eval` on the observed CFs "
          "of all quartets, which makes runs that used `propQuartets` < 1 comparable. After `bench.py rescore`, "
          "HWCD and the scores on all quartets of every version come from the same `csnaq_eval`.", ""]
    figures = []
    for sname in [s for s in suites if any(r["suite"] == s for r in recs)]:
        suite = suites[sname]
        labels = ordered_labels(recs, sname, args.baseline, chosen)
        latest = latest_runs_by_label(recs, sname)
        cases = [c["name"] for c in suite["cases"] if any(c["name"] in latest.get(l, {}) for l in labels)]
        if not labels or not cases:
            continue
        md += [f"## Suite `{sname}`", "", suite["description"], ""]
        runs_line = ", ".join(f"**{l}**: run `{next(iter(latest[l].values()))['run_id']}`" for l in labels)
        md += [f"Latest run of each version: {runs_line}.", ""]
        host = next(iter(latest[labels[-1]].values()))["host"]
        threads = next(iter(latest[labels[-1]].values()))["threads"]
        sub_host = f"{host.get('cpu', '')}, {threads} thread(s)"
        if suite["kind"] == "snaq":
            specs = [("time", "time_s", "Runtime", "seconds", True),
                     ("memory", "peak_rss_mb", "Peak memory", "MB (peak resident set size)", True),
                     ("hwcd", "hwcd", "Distance of the best network to the truth", "HWCD (lower is better)", False),
                     ("score", "score_gap", "Score of the best network on all quartets, minus the true network's",
                      "SNaQ score difference (higher is better; 0 = as good as the fitted true network)", False)]
        else:
            specs = [("fit", "fit_s", "Parameter optimization time", "seconds", True),
                     ("equations", "equations_s", "Time to build the quartet equations", "seconds", True),
                     ("memory", "peak_rss_mb", "Peak memory", "MB (peak resident set size)", True)]
        for fname, key, title, unit, log in specs:
            vals = [[metric(latest[l].get(c), key) for c in cases] for l in labels]
            svg = grouped_bar_svg(f"{title} — suite {sname}", sub_host, cases, labels, vals, unit, log=log)
            path = os.path.join(REPORTS, f"{sname}_{fname}.svg")
            if write(path, svg):
                figures.append(path)
                md += [f"![{title}]({os.path.basename(path)})", ""]
        # speedup relative to the baseline
        if args.baseline in labels and len(labels) > 1:
            base_key = "time_s" if suite["kind"] == "snaq" else "fit_s"
            others = [l for l in labels if l != args.baseline]
            vals = []
            for l in others:
                row = []
                for c in cases:
                    b, v = metric(latest[args.baseline].get(c), base_key), metric(latest[l].get(c), base_key)
                    row.append(b / v if b and v else None)
                vals.append(row)
            svg = grouped_bar_svg(f"Speedup over {args.baseline} — suite {sname}",
                                  f"{args.baseline} time ÷ version time; {sub_host}", cases, others, vals,
                                  "times faster", log=True, reference=1.0, reference_label="same speed")
            path = os.path.join(REPORTS, f"{sname}_speedup.svg")
            if write(path, svg):
                figures.append(path)
                md += [f"![Speedup]({os.path.basename(path)})", ""]
        # table
        if suite["kind"] == "snaq":
            cols = [("time (s)", "time_s"), ("peak MB", "peak_rss_mb"), ("best score", "best_score"),
                    ("median run score", "median_run_score"), ("score on all quartets", "score_all_quartets"),
                    ("truth score on all quartets", "score_truth_all_quartets"), ("HWCD", "hwcd"),
                    ("hybrids", "nhybrids")]
        else:
            cols = [("fit (s)", "fit_s"), ("equations (s)", "equations_s"), ("peak MB", "peak_rss_mb"),
                    ("score", "score"), ("quartets", "nquartets"), ("parameters", "nparams")]
        md += ["| case | version | " + " | ".join(c for c, _ in cols) + " |",
               "|---|---|" + "---|" * len(cols)]
        for c in cases:
            for l in labels:
                rec = latest[l].get(c)
                if rec is None:
                    continue
                if rec["metrics"].get("failed"):
                    md.append(f"| {c} | {l} | failed |" + " |" * (len(cols) - 1))
                    continue
                cells = []
                for _, k in cols:
                    v = metric(rec, k)
                    cells.append("" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v)))
                md.append(f"| {c} | {l} | " + " | ".join(cells) + " |")
        md.append("")
        # history: geometric-mean speedup over the baseline, for every run of the suite
        if any(r["label"] == args.baseline and r["suite"] == sname for r in recs):
            base_key = "time_s" if suite["kind"] == "snaq" else "fit_s"
            base = latest[args.baseline]
            runs = {}
            for r in sorted(recs, key=lambda r: r["run_id"]):
                if r["suite"] == sname and r["label"] != args.baseline:
                    runs.setdefault(r["run_id"], []).append(r)
            if runs:
                xl, sp, hw = [], [], []
                for rid, rs in runs.items():
                    ratios = [metric(base.get(r["case"]), base_key) / metric(r, base_key) for r in rs
                              if metric(r, base_key) and metric(base.get(r["case"]), base_key)]
                    xl.append(f"{rs[0]['label']} ({rid[:8]})")
                    sp.append(geomean(ratios))
                svg = line_svg(f"Speedup over {args.baseline} across versions — suite {sname}",
                               "geometric mean over the suite's cases, for every recorded run (oldest first)",
                               xl, ["geometric-mean speedup"], [sp], "times faster", log=True, reference=1.0,
                               reference_label=f"{args.baseline} speed")
                path = os.path.join(REPORTS, f"{sname}_history.svg")
                if write(path, svg):
                    figures.append(path)
                    md += [f"![History]({os.path.basename(path)})", ""]
    with open(os.path.join(REPORTS, "REPORT.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"wrote {os.path.relpath(os.path.join(REPORTS, 'REPORT.md'), ROOT)} and {len(figures)} figure(s)")


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, default_suite):
        sp.add_argument("--suite", default=default_suite, help=f"suite in benchmarks/suites.json (default: {default_suite})")
        sp.add_argument("--cases", help="comma-separated subset of the suite's cases")
        sp.add_argument("--label", help="name of this version in the history (default: implementation@commit)")
        sp.add_argument("--threads", type=int, default=1, help="threads (default: 1)")
        sp.add_argument("--notes", help="free-text notes stored with the results")
        sp.add_argument("--eval-binary", default=DEFAULT_EVAL, help="csnaq_eval, for accuracy metrics")
        sp.add_argument("--repeats", type=int, default=3, help="repeats of each fit case (minimum time is kept)")

    r = sub.add_parser("run", help="benchmark a snaq binary (C++)")
    common(r, "standard")
    r.add_argument("--binary", default=DEFAULT_SNAQ, help="snaq binary to benchmark (default: build/snaq)")
    r.add_argument("--fit-binary", default=DEFAULT_FIT, help="csnaq_bench_fit (for the fit suite)")
    r.add_argument("--extra-args", help='extra snaq arguments, e.g. "--my-experimental-flag 3"')

    j = sub.add_parser("run-julia", help="benchmark SNaQ.jl (the baseline)")
    common(j, "standard")
    j.add_argument("--julia", default="julia", help="julia executable")

    rp = sub.add_parser("report", help="draw figures and write benchmarks/reports/REPORT.md")
    rp.add_argument("--baseline", default="SNaQ.jl", help="label of the reference version (default: SNaQ.jl)")
    rp.add_argument("--labels", help="comma-separated labels to compare (default: baseline + 3 latest)")

    sub.add_parser("list", help="list the recorded benchmark runs")

    rs = sub.add_parser("rescore", help="recompute the accuracy metrics of all recorded snaq runs with one csnaq_eval")
    rs.add_argument("--eval-binary", default=DEFAULT_EVAL, help="csnaq_eval to use (default: build/csnaq_eval)")
    rs.add_argument("--threads", type=int, default=1, help="threads (default: 1)")

    args = p.parse_args()
    if args.cmd == "run":
        run_suite(args, "cpp")
    elif args.cmd == "run-julia":
        args.binary = args.fit_binary = None
        args.extra_args = None
        if args.label is None:
            args.label = "SNaQ.jl"
        run_suite(args, "julia")
    elif args.cmd == "report":
        report(args)
    elif args.cmd == "list":
        list_runs(args)
    elif args.cmd == "rescore":
        rescore(args)


if __name__ == "__main__":
    main()
