import sys
import subprocess
import toml

TOML_PATH = "pyproject.toml"

def get_required_versions():
    data = toml.load(TOML_PATH)
    python_req = data.get("project", {}).get("requires-python", None)
    deps = data.get("project", {}).get("dependencies", [])
    required = {}
    if python_req:
        required["python"] = python_req
    for dep in deps:
        if ">=" in dep or "==" in dep:
            pkg, ver = dep.replace(" ", "").split(">=" if ">=" in dep else "==")
            required[pkg.lower()] = ver
        else:
            required[dep.lower()] = None
    return required

def get_installed_versions(required):
    installed = {}
    py_ver = ".".join(map(str, sys.version_info[:3]))
    installed["python"] = py_ver
    for pkg in required:
        if pkg == "python":
            continue
        try:
            out = subprocess.check_output([sys.executable, "-m", "pip", "show", pkg], text=True)
            for line in out.splitlines():
                if line.startswith("Version:"):
                    installed[pkg] = line.split(":", 1)[1].strip()
        except Exception:
            installed[pkg] = None
    return installed

def compare_versions(required, installed):
    mismatches = []
    for pkg, req_ver in required.items():
        inst_ver = installed.get(pkg)
        if req_ver and inst_ver and not inst_ver.startswith(req_ver):
            mismatches.append(f"{pkg}: required {req_ver}, found {inst_ver}")
        elif req_ver and not inst_ver:
            mismatches.append(f"{pkg}: required {req_ver}, not installed")
    return mismatches

def main():
    required = get_required_versions()
    installed = get_installed_versions(required)
    mismatches = compare_versions(required, installed)
    if mismatches:
        print("Version mismatches found:")
        for m in mismatches:
            print(m)
        sys.exit(1)
    print("All versions match.")
    sys.exit(0)

if __name__ == "__main__":
    main()