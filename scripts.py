import os
import re
import pexpect
import pandas as pd
import io
from tabulate import tabulate

PROG_PATH = "./build/main_snap"
ORDERS_WITH_TIMEOUT = [3, 4, 5, 6, 7]

def find_limited_files_per_topdir(root_dir, limit):
    """
    For each immediate subdirectory in root_dir, find the first file
    in any of its nested subdirectories (recursively), and stop after
    reaching the specified limit of files.
    
    :param root_dir: The root directory to scan
    :param limit: Maximum number of files to find
    """
    root_dir = os.path.abspath(os.path.expanduser(root_dir))
    count = 0
    res = []
    
    # Iterate over immediate subdirectories of root_dir
    for entry in os.scandir(root_dir):
        if entry.is_dir():
            if count >= limit:
                break
            big_dir = entry.path
            found = False
            # Walk recursively but stop after finding the first file
            for subdir, _, files in os.walk(big_dir):
                if files:
                    first_file = os.path.join(subdir, files[0])
                    print(f"File {count+1}: {first_file}")
                    res.append(first_file)
                    found = True
                    count += 1
                    break  # stop at first file
            if not found:
                print(f"No files found in '{big_dir}'")
    return res

def find_all_files_per_topdir(root_dir, limit=None):
    """
    Collect all files under root_dir (including those directly inside it)
    and under each of its subdirectories, recursively.
    Optionally stop after reaching 'limit' files total.

    :param root_dir: The root directory to scan
    :param limit: Maximum number of files to find (None = no limit)
    """
    root_dir = os.path.abspath(os.path.expanduser(root_dir))
    count = 0
    res = []

    # First: include files directly under root_dir
    for entry in os.scandir(root_dir):
        if entry.is_file():
            if limit is not None and count >= limit:
                return res
            res.append(entry.path)
            count += 1

    # Then: recurse into immediate subdirectories
    for entry in os.scandir(root_dir):
        if entry.is_dir():
            for subdir, _, files in os.walk(entry.path):
                for fname in files:
                    if limit is not None and count >= limit:
                        return res
                    file_path = os.path.join(subdir, fname)
                    res.append(file_path)
                    count += 1

    return res

# Define CSV schema (all columns, so new files get them initialized)
COLUMNS = [
    "file", "rows", "cols", "nnz", "initial_diag",
    "rcm_final_diag", "rcm_time",
    "gorder_w10_final_diag", "gorder_w10_time",
    "gorder_w50_final_diag", "gorder_w50_time",
    "gorder_w100_final_diag", "gorder_w100_time",
    "hs_final_diag", "hs_time",
    "rcmhs_final_diag", "rcmhs_time",
    "gorderhs_w10_final_diag", "gorderhs_w10_time",
    "gorderhs_w50_final_diag", "gorderhs_w50_time",
    "gorderhs_w100_final_diag", "gorderhs_w100_time",
    "hsbf_final_diag", "hsbf_time"
]

def update_csv(csv_file, abs_file, order_id, window_size, parsed):
    """
    Update the CSV file with results for (abs_file, order_id, window_size).
    If the row doesn't exist, create it.
    """
    # Load or initialize dataframe
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
    else:
        df = pd.DataFrame(columns=COLUMNS)

    # Ensure all expected columns exist
    for col in COLUMNS:
        if col not in df.columns:
            df[col] = None

    # Find row index for this file
    if abs_file in df["file"].values:
        idx = df.index[df["file"] == abs_file][0]
    else:
        # New row
        idx = len(df)
        df.loc[idx] = {col: None for col in COLUMNS}
        df.loc[idx, "file"] = abs_file
        df.loc[idx, "rows"] = parsed.get("rows")
        df.loc[idx, "cols"] = parsed.get("cols")
        df.loc[idx, "nnz"] = parsed.get("nnz")
        df.loc[idx, "initial_diag"] = parsed.get("initial_diagonal")

    df[COLUMNS] = df[COLUMNS].fillna("").astype(str)

    # Map order_id + window_size → column names
    if order_id == 1:
        df.loc[idx, "rcm_final_diag"] = parsed.get("final_diagonal")
        df.loc[idx, "rcm_time"] = parsed.get("time_rcm")

    elif order_id == 2:
        if window_size == 100:
            df.loc[idx, "gorder_w100_final_diag"] = parsed.get("final_diagonal")
            df.loc[idx, "gorder_w100_time"] = parsed.get("time_gorder")
        elif window_size == 10:
            df.loc[idx, "gorder_w10_final_diag"] = parsed.get("final_diagonal")
            df.loc[idx, "gorder_w10_time"] = parsed.get("time_gorder")
        elif window_size == 50:
            df.loc[idx, "gorder_w50_final_diag"] = parsed.get("final_diagonal")
            df.loc[idx, "gorder_w50_time"] = parsed.get("time_gorder")

    elif order_id == 3:
        df.loc[idx, "hs_final_diag"] = parsed.get("final_diagonal")
        df.loc[idx, "hs_time"] = parsed.get("time_hs")

    elif order_id == 4:
        df.loc[idx, "rcmhs_final_diag"] = parsed.get("final_diagonal")
        df.loc[idx, "rcmhs_time"] = parsed.get("time_rcm") + parsed.get("time_hs")

    elif order_id == 5:
        if window_size == 100:
            df.loc[idx, "gorderhs_w100_final_diag"] = parsed.get("final_diagonal")
            df.loc[idx, "gorderhs_w100_time"] = parsed.get("time_gorder") + parsed.get("time_hs")
        elif window_size == 10:
            df.loc[idx, "gorderhs_w10_final_diag"] = parsed.get("final_diagonal")
            df.loc[idx, "gorderhs_w10_time"] = parsed.get("time_gorder") + parsed.get("time_hs")
        elif window_size == 50:
            df.loc[idx, "gorderhs_w50_final_diag"] = parsed.get("final_diagonal")
            df.loc[idx, "gorderhs_w50_time"] = parsed.get("time_gorder") + parsed.get("time_hs")
    
    elif order_id == 7:
        df.loc[idx, "hsbf_final_diag"] = parsed.get("final_diagonal")
        df.loc[idx, "hsbf_time"] = parsed.get("time_hsbf")


    # Save CSV back
    df.to_csv(csv_file, index=False)

def run_with_timeout(cmd, timeout=600):
    child = pexpect.spawn(" ".join(cmd), encoding="utf-8", timeout=None)

    # capture everything the child prints
    logbuf = io.StringIO()
    child.logfile_read = logbuf

    try:
        child.expect("Read Matrix", timeout=None)
    except pexpect.TIMEOUT:
        print(f"{' '.join(cmd)} did not print 'Read Matrix'")
        child.terminate(force=True)
        return logbuf.getvalue()

    try:
        child.expect(pexpect.EOF, timeout=timeout)
    except pexpect.TIMEOUT:
        print(f"Timeout expired for {' '.join(cmd)}, sending key press...")
        child.sendline("")
        try:
            child.expect("Result Matrix")
            child.expect(pexpect.EOF)
        except pexpect.TIMEOUT:
            child.terminate(force=True)

    stdout = logbuf.getvalue()
    return stdout



def run_without_timeout(cmd):
    """
    Run command without timeout, wait until it finishes.
    Returns stdout as a string.
    """
    child = pexpect.spawn(" ".join(cmd), encoding="utf-8")

    # Wait until program finishes
    child.expect(pexpect.EOF)

    stdout = child.before
    return stdout

def run_order_id(file, order_id, window_size, timeout=600):
    if order_id in [2, 5]:
        cmd = [PROG_PATH, file, str(order_id), str(window_size)]
        print(f"Running: {' '.join(cmd)}")
        if order_id == 5:
            stdout = run_with_timeout(cmd, timeout=timeout)
        else:
            stdout = run_without_timeout(cmd)
    else:
        cmd = [PROG_PATH, file, str(order_id)]
        print(f"Running: {' '.join(cmd)}")
        if order_id in ORDERS_WITH_TIMEOUT:
            stdout = run_with_timeout(cmd, timeout=timeout)
        else:
            stdout = run_without_timeout(cmd)
    return stdout


def parse_output(stdout, order_id):
    """
    Parse program stdout and extract diagonal counts and timings.
    Returns dict with parsed data.
    """
    initial_diag = None
    final_diag = None
    intermediate_diag = None

    time_rcm = None
    time_gorder = None
    time_hs = None
    time_hsbf = None

    for line in stdout.splitlines():
        line = line.strip()

        if line.startswith("Initial diagonal count:"):
            try:
                initial_diag = int(line.split(":")[1].strip())
            except:
                pass

        elif re.match(r"^Diagonal count:\s*\d+", line):  
            # Final diagonal (not initial or intermediate)
            try:
                final_diag = int(line.split(":")[1].strip())
            except:
                pass

        elif line.startswith("Intermediate diagonal count after RCM:"):
            try:
                intermediate_diag = int(line.split(":")[1].strip())
            except:
                pass

        elif line.startswith("Intermediate diagonal count after GOrder:"):
            try:
                intermediate_diag = int(line.split(":")[1].strip())
            except:
                pass

        # timings
        m = re.search(r"RCM Ordering took (\d+) ms", line)
        if m:
            time_rcm = int(m.group(1))

        m = re.search(r"GOrder took (\d+) ms", line)
        if m:
            time_gorder = int(m.group(1))

        m = re.search(r"HS Ordering took (\d+) ms", line)
        if m:
            time_hs = int(m.group(1))

        m = re.search(r"HSBF Ordering took (\d+) ms", line)
        if m:
            time_hsbf = int(m.group(1))

    return {
        "initial_diagonal": initial_diag,
        "intermediate_diagonal": intermediate_diag,
        "final_diagonal": final_diag,
        "time_rcm": time_rcm,
        "time_gorder": time_gorder,
        "time_hs": time_hs,
        "time_hsbf": time_hsbf
    }

def run_main_on_file(file, csv_file, interrupt_timeout=600):
    abs_file = os.path.abspath(os.path.expanduser(file))
    rows, cols, nnz = None, None, None

    with open(abs_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue  # skip comments
            # first non-comment line after banner is the size info
            parts = line.split()
            if len(parts) == 3:
                rows, cols, nnz = map(int, parts)
                break
            else:
                raise ValueError("Unexpected format in .mtx file header")

    order_ids = [1, 2, 3, 4, 5, 7]
    window_sizes = [10, 50, 100]

    for order_id in order_ids:
        if order_id in [2, 5]:
            for w in window_sizes:
                stdout = run_order_id(abs_file, order_id, w, timeout_duration)
                parsed = parse_output(stdout, order_id)

                parsed = {
                    "rows": rows,
                    "cols": cols,
                    "nnz": nnz,
                    **parsed
                }
                
                update_csv(csv_file, os.path.splitext(os.path.basename(abs_file.split()[-1]))[0], order_id, w, parsed)

        else:
            stdout = run_order_id(abs_file, order_id, 0, timeout_duration)
            parsed = parse_output(stdout, order_id)

            parsed = {
                "rows": rows,
                "cols": cols,
                "nnz": nnz,
                **parsed
            }

            update_csv(csv_file, os.path.splitext(os.path.basename(abs_file.split()[-1]))[0], order_id, 0, parsed)

def run_main_on_files(files, csv_file, interrupt_timeout=600):
    for file in files:
        run_main_on_file(file, csv_file, interrupt_timeout)
    
def run_one_ordering(file, csv_file, order_id, w, timeout):
    abs_file = os.path.abspath(os.path.expanduser(file))
    rows, cols, nnz = None, None, None

    with open(abs_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue  # skip comments
            # first non-comment line after banner is the size info
            parts = line.split()
            if len(parts) == 3:
                rows, cols, nnz = map(int, parts)
                break
            else:
                raise ValueError("Unexpected format in .mtx file header")
    
    stdout = run_order_id(file, order_id, w, timeout)
    parsed = parse_output(stdout, order_id)

    parsed = {
        "rows": rows,
        "cols": cols,
        "nnz": nnz,
        **parsed
    }
    
    update_csv(csv_file, os.path.splitext(os.path.basename(abs_file.split()[-1]))[0], order_id, w, parsed)

def write_table_to_file(filename, csv_file):
    # Read your CSV
    df = pd.read_csv(csv_file)

    # Convert dataframe to table string
    table_str = tabulate(df, headers="keys", tablefmt="grid")

    # Save to a text file
    with open(filename, "w") as f:
        f.write(table_str)

    print(f"Table written to {filename}")

def run_one_ordering_on_dir(files, csv_file, order_id, w, interrupt_timeout=600):
    for file in files:
        run_one_ordering(file, csv_file, order_id, w, interrupt_timeout)

if __name__ == "__main__":
    print("Choose action: ")
    print("    1- Write table from csv file")
    print("    2- Run program on a matrix")
    print("    3- Run program on matrices in dir")
    print("    4- Run one ordering on a matrix")
    print("    5- Run one ordering on matrices in a dir")

    action_choice = int(input("Choice: "))

    if action_choice == 1:
        output_filename = input("Enter output filename: ")
        input_filename = input("Enter csv input filename: ")
        write_table_to_file(output_filename, input_filename)
    elif action_choice == 2:
        mtx = input("Enter path to matrix: ")
        csv_file = input("Enter csv filename: ")
        timeout_duration = int(input("Enter timeout for HSOrder in sec: "))
        run_main_on_file(mtx, csv_file, timeout_duration)
    elif action_choice == 3:
        dirname = input("Enter path to dir: ")
        csv_file = input("Enter path to csv file: ")
        timeout_duration = int(input("Enter timeout for HSOrder in sec: "))
        files = find_all_files_per_topdir(dirname)
        run_main_on_files(files, csv_file, timeout_duration)
    elif action_choice == 4:
        mtx = input("Enter path to matrix: ")
        csv_file = input("Enter csv filename: ")
        order_id = int(input("Enter order id: "))
        if order_id in [2, 5]:
            w = int(input("Enter window size (10, 50, 100): "))
            assert w in [10, 50, 100]
        else:
            w = 0
        if order_id in ORDERS_WITH_TIMEOUT:
            timeout_duration = int(input("Enter timeout for HSOrder in sec: "))
        else:
            timeout_duration = 0
        run_one_ordering(mtx, csv_file, order_id, w, timeout_duration)
    elif action_choice == 5:
        dirname = input("Enter path to dir: ")
        csv_file = input("Enter path to csv file: ")
        order_id = int(input("Enter order id: "))
        if order_id in [2, 5]:
            w = int(input("Enter window size (10, 50, 100): "))
            assert w in [10, 50, 100]
        else:
            w = 0
        if order_id in ORDERS_WITH_TIMEOUT:
            timeout_duration = int(input("Enter timeout for HSOrder in sec: "))
        else:
            timeout_duration = 0
        files = find_all_files_per_topdir(dirname)
        run_one_ordering_on_dir(files, csv_file, order_id, w, timeout_duration)
