from __future__ import annotations
import hashlib, shlex, subprocess, sys, time
from pathlib import Path

def eprint(*args, **kwargs): print(*args, file=sys.stderr, flush=True, **kwargs)
def cmd_text(cmd): return " ".join(shlex.quote(str(x)) for x in cmd)
def run(cmd, *, stdout=None, stderr=None, check=True, text=True, cwd=None):
    eprint('+', cmd_text(cmd))
    return subprocess.run([str(x) for x in cmd], stdout=stdout, stderr=stderr, check=check, text=text, cwd=cwd)
def capture(cmd, *, cwd=None):
    return run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd).stdout
def retry_run(cmd, tries=4, sleep_s=4, **kwargs):
    """Retry a subprocess command with increasing delay; raise after the final failure."""
    last=None
    for attempt in range(1,tries+1):
        eprint(f"+ [attempt {attempt}/{tries}] {cmd_text(cmd)}")
        p=subprocess.run([str(x) for x in cmd], **kwargs)
        if p.returncode==0: return p
        last=p
        if attempt<tries: time.sleep(min(sleep_s*attempt,20))
    raise RuntimeError(f"Command failed after {tries} attempts (exit={last.returncode if last else 'NA'}): {cmd_text(cmd)}")

def retry_shell_capture(command: str, tries=4, sleep_s=4):
    last=None
    for attempt in range(1, tries+1):
        eprint(f'+ [attempt {attempt}/{tries}] {command}')
        p=subprocess.run(['bash','-lc',command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.returncode==0 and p.stdout.strip(): return p.stdout
        last=p
        if attempt<tries: time.sleep(sleep_s*attempt)
    msg=(last.stderr if last else '')[-3000:]
    raise RuntimeError(f'Command failed after {tries} attempts:\n{command}\n{msg}')
def ensure_executable(name):
    from shutil import which
    if which(name) is None: raise RuntimeError(f'Required executable not found in PATH: {name}')
def md5_file(path: Path, chunk=8*1024*1024):
    h=hashlib.md5()
    with path.open('rb') as f:
        while True:
            b=f.read(chunk)
            if not b: break
            h.update(b)
    return h.hexdigest()
def now_iso(): return time.strftime('%Y-%m-%dT%H:%M:%S%z')
