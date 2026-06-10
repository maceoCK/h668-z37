"""Modal: short annealing runs with best-state residual dumps (floor fingerprint analysis)."""
import modal
app = modal.App("h668-sc-residual")
image = (modal.Image.debian_slim(python_version="3.11").apt_install("g++")
         .add_local_file("sc_anneal.cpp", "/root/sc_anneal.cpp")
         .add_local_file("sc_tpl_221_0.txt", "/root/sc_tpl_221_0.txt")
         .add_local_file("sc_tpl_221_1.txt", "/root/sc_tpl_221_1.txt")
         .add_local_file("sc_tpl_422_2.txt", "/root/sc_tpl_422_2.txt")
         .add_local_file("sc_tpl_422_3.txt", "/root/sc_tpl_422_3.txt")
         .add_local_file("sc_tpl_422_4.txt", "/root/sc_tpl_422_4.txt")
         .add_local_file("sc_tpl_422_5.txt", "/root/sc_tpl_422_5.txt"))

@app.function(image=image, cpu=1.0, timeout=3000, retries=2)
def run_one(tpl: str, seed: int, time_s: int):
    import subprocess
    subprocess.run(["g++", "-O3", "-o", "/root/sc", "/root/sc_anneal.cpp"], check=True)
    pr = subprocess.run(["/root/sc", f"/root/{tpl}", str(seed), str(time_s)],
                        capture_output=True, text=True, timeout=time_s + 300)
    keep = [ln for ln in pr.stdout.splitlines()
            if ln.startswith(("BESTSTATE", "REPBITS", "RES ", "RESULT", "SOLUTION", "D "))]
    return {"tpl": tpl, "seed": seed, "dump": "\n".join(keep)}

@app.local_entrypoint()
def main(time_s: int = 1800):
    import json
    tpls = ["sc_tpl_221_0.txt", "sc_tpl_221_1.txt", "sc_tpl_422_2.txt",
            "sc_tpl_422_3.txt", "sc_tpl_422_4.txt", "sc_tpl_422_5.txt"]
    jobs = [(t, s, time_s) for t in tpls for s in (101, 102)]
    out = []
    for r in run_one.starmap(jobs):
        out.append(r)
        first = r["dump"].splitlines()[0] if r["dump"] else "?"
        print(f"  {r['tpl']} seed {r['seed']}: {first}", flush=True)
    with open("sc_residuals.json", "w") as f:
        json.dump(out, f, indent=1)
    print("=== RESIDUAL DUMPS SAVED ===")
