#!/usr/bin/env bash
# Pre-flight check: is this computer ready to run the pipeline? Prints what is wrong and how to fix it.
ok=0; bad=0
pass() { echo "  [ok]   $*"; ok=$((ok+1)); }
fail() { echo "  [FAIL] $*"; bad=$((bad+1)); }
warn() { echo "  [warn] $*"; }
echo "Checking this computer..."
if command -v java >/dev/null 2>&1; then
    jv=$(java -version 2>&1 | head -1 | sed -E 's/.*"([0-9]+)[^"]*".*/\1/')
    [ "${jv:-0}" -ge 17 ] 2>/dev/null && pass "java $jv" || fail "java is $jv, need 17 or newer:  sudo apt install openjdk-17-jre-headless"
else fail "java not found:  sudo apt install openjdk-17-jre-headless"; fi
if command -v nextflow >/dev/null 2>&1; then pass "nextflow $(nextflow -version 2>/dev/null | grep -oE 'version [0-9.]+' | head -1)"
else fail "nextflow not found:  curl -s https://get.nextflow.io | bash && mkdir -p ~/.local/bin && mv nextflow ~/.local/bin/"; fi
if command -v docker >/dev/null 2>&1; then
    if docker info >/dev/null 2>&1; then pass "docker daemon is running"
    else fail "docker is installed but not running or you lack permission. WSL: start Docker Desktop on Windows (and enable WSL integration), or: sudo service docker start; sudo usermod -aG docker \$USER (then log out/in)"; fi
else fail "docker not found: install Docker Desktop (Windows/Mac) or docker.io (Linux)"; fi
command -v python3 >/dev/null 2>&1 && pass "python3" || fail "python3 not found:  sudo apt install python3"
mem=$(awk '/MemTotal/{printf "%d", $2/1024/1024}' /proc/meminfo 2>/dev/null)
if [ -n "$mem" ]; then [ "$mem" -ge 12 ] && pass "memory ${mem} GB" || warn "memory ${mem} GB: assemblies may be killed; run with --max_memory '${mem}.GB' and expect retries"; fi
cpus=$(nproc 2>/dev/null); [ -n "$cpus" ] && pass "cpus $cpus"
free=$(df -BG "${1:-.}" 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}')
if [ -n "$free" ]; then [ "$free" -ge 20 ] && pass "disk ${free} GB free" || warn "only ${free} GB free: Docker images need ~6 GB, results ~1-3 GB per sample"; fi
case "$(readlink -f "${1:-.}")" in /mnt/[a-z]/*) warn "input is on the Windows drive (/mnt/...): reading is slow. Copy the reads to a folder under /home first.";; esac
if grep -qi microsoft /proc/version 2>/dev/null; then pass "WSL2 detected"; fi
echo
if [ "$bad" -gt 0 ]; then echo "$bad problem(s) found. Fix them and run again."; exit 1; fi
echo "All checks passed."
