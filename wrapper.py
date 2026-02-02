"""
Luke Cirne
ColabFold Wrapper
Template Cycling with Experimental Distance Restraint Data
Ma Lab
"""
import os
import sys
import subprocess
import shutil
import json
import getpass
import numpy as np
import data_engine as engine

def initialize_project(jobs):
    """
    Initializes a ColabFold project by gathering user input,
    creating necessary variables, generating a shell script,
    and appending metadata to a JSON log.

    Args:
        jobs (str): Path to the JSON file tracking job metadata.

    Returns:
        str: Path to the generated shell script for running ColabFold.
    """
    key_values = []
    # Obtain username
    try:
        username = getpass.getuser()
        print(f">>> LOGGING AS: {username}")
    except OSError:
        print(">>> COULD NOT OBTAIN USERNAME")
        print(">>> LOGGING AS: DEFAULT")
        username = "DEFAULT"
    key_values.append(("user", username))

    # Obtain job ID
    try:
        with open(jobs, "r") as file:
            print(">>> READING JSON")
            job_dict = json.load(file)
            user_JIDs = [int(job["JID"]) for job in job_dict["jobs"] if job["user"] == username]
            last_user_JID = max(user_JIDs) if user_JIDs else 0
            current_JID = last_user_JID + 1
            key_values.append(("JID", current_JID))
    except FileNotFoundError:
        current_JID = 0
        key_values.append(("JID", current_JID))

    # Obtain input file
    while True:
        input_file = input("Input file name in the format 'name.fasta': ")
        if os.path.isfile(f"./{input_file}"):
            break
        else:
            print("###### INVALID FILEPATH ######")
    key_values.append(("input_file", input_file))

    # Obtain template directory
    while True:
        temp_dir = input("Input template directory: ")
        if os.path.isdir(temp_dir):
            break
        else:
            print("###### INVALID DIRECTORY ######")
    key_values.append(("temp_dir", temp_dir))

    # Obtain num recycles
    while True:
        num_c = input("Desired number of recycles (integer) (max 5, min 0): ")
        try:
            if 0 <= int(num_c) < 6:
                break
            else:
                print("###### Invalid input ######")
        except ValueError:
                print("###### Invalid input ######")
    key_values.append(("num_recycles", num_c))    

    # Obtain num seeds
    while True:
        num_s = input("Desired number of seeds (integer) (min 1): ")
        try:
            if 0 < int(num_s):
                num_s = int(num_s)
                break
            else:
                print("###### Invalid input #######")
        except ValueError:
                print("###### Invalid input #######")
    key_values.append(("num_s", num_s))

    # Obtain value for n
    while True:
        n = input("Desired value for n (integer) (min 10): ")
        try:
            if 9 < int(n):
                n = int(n)
                break
            else:
                print("###### Invalid input #######")
        except ValueError:
                print("###### Invalid input #######")
    key_values.append(("n", n))

    # Obtain value for m_e_msa
    while True:
        m_e_msa = input("Desired value for max number of msa's (integer) (min 1): ")
        try:
            if 1 < int(m_e_msa):
                m_e_msa = int(m_e_msa)
                break
            else:
                print("###### Invalid input #######")
        except ValueError:
                print("###### Invalid input #######")
    key_values.append(("m_e_msa", m_e_msa))

    m_msa = m_e_msa // 2
    key_values.append(("m_msa", m_msa))

    # Variable for full output directory by user, job id, and max msa's
    outputdir = f"{username}{current_JID}mm{m_msa}"
    key_values.append(("outputdir", outputdir))

    # Create shell script to run colabfold_batch
    current_dir = os.path.dirname(os.path.abspath(__file__))
    script_name = "wrapper.sh"
    script_path = os.path.join(current_dir, script_name)
    # Writing script
    script_content = f"""#!/bin/bash
JID={current_JID}
num_c={num_c}
seed=1
num_s={num_s}
m_e_msa={m_e_msa}
m_msa={m_msa}
inputfile=./{input_file}
outputdir={outputdir}
temp_dir={temp_dir}

colabfold_batch --pair-mode unpaired_paired --templates \\
--msa-mode mmseqs2_uniref_env \\
--custom-template-path $temp_dir \\
--max-msa $m_msa:$m_e_msa \\
--use-dropout \\
--num-seeds $num_s \\
--num-recycle $num_c \\
$inputfile $outputdir
    """
    # Create shell script to execute ColabFold
    with open(script_path, 'w') as file:
        print(">>> WRITING SHELL SCRIPT")
        file.write(script_content)
    append_jobs_json("jobs.json", key_values)
    return script_path


def delete_directory(dir_path):
    """
    Deletes a directory and all its contents.

    Args:
        dir_path (str): The path to the directory to delete.

    Returns:
        bool: True if the directory was successfully deleted, False otherwise.
    """
    try:
        shutil.rmtree(dir_path)
        return True
    except OSError:
        return False        


def clear_directory(dir_path):
    """
    Clears the contents of a directory without deleting the directory itself.

    Args:
        dir_path (str): Path to the directory to be cleared.
    """
    for item in os.listdir(dir_path):
        item_path = os.path.join(dir_path, item)
        if os.path.isfile(item_path):
            os.remove(item_path)  # Remove files
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)


def append_jobs_json(jobs, key_values):
    """
    Appends job metadata to a JSON file. If the file does not exist,
    it creates a new one.

    Args:
        jobs (str): Path to the JSON file.
        key_values (list): List of (key, value) tuples to append.
    """
    # Convert list of tuples into a proper dictionary
    new_entry = {key: value for key, value in key_values}
    try:
        with open(jobs, "r") as f:
            job_dict = json.load(f)
    except FileNotFoundError:
        # Create new JSON structure if file doesn't exist
        job_dict = {"jobs": []}

    # Append the new dictionary entry
    job_dict["jobs"].append(new_entry)

    with open(jobs, "w") as f:
        print(f">>> APPENDING JSON TO {jobs}")
        json.dump(job_dict, f, indent=4)
        print("###### COMPLETE ######")


def append_mods_json(mods_file, mods_dict):
    """
    Append distribution modifications from the distribution building
    process to the json logs.

    Args:
        mods_file (str): Path to JSON file.
        mods_json (dict): Dictionary containing json to append.
    """
    try:
        with open(mods_file, "r") as f:
            mods_json = json.load(f)
    except FileNotFoundError as e:
        print(f">>> EXCEPTION WHEN APPENDING TO {mods_file}: {e}")
    except json.JSONDecodeError:
        mods_json = {"mods": []}
    mods_json["mods"].append(mods_dict)

    with open(mods_file, "w") as f:
        print(f">>> APPENDING JSON TO {mods_file}")
        json.dump(mods_json, f, indent=4)
        print("###### COMPLETE ######")


def get_from_current_job(jobs_file, items) -> list:
    try:
        with open(jobs_file, "r") as file:
            jobs_dict = json.load(file)
            current_job_info = list(jobs_dict["jobs"][-1].values())
            requested_items = []
            for item in items:
                match item:
                    case "user":
                        requested_items.append(current_job_info[0])
                    case "JID":
                        requested_items.append(current_job_info[1])
                    case "input_file":
                        requested_items.append(current_job_info[2])
                    case "temp_dir":
                        requested_items.append(current_job_info[3])
                    case "num_recycles":
                        requested_items.append(current_job_info[4])
                    case "num_s":
                        requested_items.append(current_job_info[5])
                    case "n":
                        requested_items.append(current_job_info[6])
                    case "m_e_msa":
                        requested_items.append(current_job_info[7])
                    case "m_msa":
                        requested_items.append(current_job_info[8])
                    case "outputdir":
                        requested_items.append(current_job_info[9])
        return requested_items
    except FileNotFoundError:
        print("###### JSON FILE NOT FOUND ######")
    return 0


def filter_output(run_number, jobs, script_path, n):
    """
    Filters PDB output files based on proximity to target distances.
    Updates template directory for next ColabFold iteration accordingly.

    Args:
        run_number (int): The current iteration number.
        jobs (str): Path to the job metadata JSON file.
        script_path (str): Path to the ColabFold execution script.
    """
    # Load json and obtain outputdir and temp_dir
    outputdir, temp_dir = get_from_current_job(jobs, ["outputdir", "temp_dir"])

    current_dir = os.getcwd()
    #colabfold_output = os.listdir(f"{current_dir}/{outputdir}")
    colabfold_output = os.listdir(f"{outputdir}")

   # Loop through files and run DistanceFinder.py on each
    distances = {}
    for file in colabfold_output:
        if file.endswith(".pdb"):
            distances[file] = float(run_distance_finder(f"{outputdir}/{file}", "100", "473"))

    distances_to_convert = np.array(list(distances.values()))
    e_conversions = engine.compute_E(distances_to_convert)
    i = 0
    for filename, distance in distances.items():
        distances[filename] = e_conversions[i]
        i += 1

    # Finished implementation:
    # write algorithm to determine which files to extract from distances
    # and add to included_distances to fit a normal distribution
    # remember to normalize data points
    # duplicate templates if necessary to fit proper distribution.
    y_exp = 0.291
    sigma = 0.083
    included_distances, bins, bin_centers, mod_count = engine.build_distribution(file_eff_dict=distances, mean=y_exp, std=sigma, n=n)

    # Save original distances using bins from build_distribution
    plot_and_save_distances(distances, run_number, bin_centers, n)
    # If included_distances dictionary is still empty after checks,
    # proceed to next iteration with user provided templates 
    if not included_distances:
        print("###### NO VALID TEMPLATES PRODUCED ######")
        update_temp_dir(script_path, temp_dir)
    else:
        temp_dir = f"iteration{run_number + 1}"
        try:
            os.mkdir("iterations")
            print(">>> CREATING DIRECTORY FOR BEST TEMPLATES")
        except FileExistsError:
            print(">>> APPENDING TO ITERATIONS DIRECTORY")

        try:
            os.mkdir(f"iterations/{temp_dir}")
        except FileExistsError:
            shutil.rmtree(f"iterations/{temp_dir}")
            os.mkdir(f"iterations/{temp_dir}")

        template_number = 0
        for filename, distance in included_distances.items():
            # Logic to check if a filename is duplicated or not
            # if so, cp the original file with new name and add to temp_dir
            # if not, just cp original file to temp_dir
            if "_dupe" in filename:
                seperator = "_dupe"
                filename_parts = filename.split(seperator, 1)
                original_filename = f"{filename_parts[0]}.pdb"
                subprocess.run(["cp", f"{outputdir}/{original_filename}", f"{outputdir}/{filename}"])

            subprocess.run(["cp", f"{outputdir}/{filename}", f"iterations/{temp_dir}"])
            template_number_str = f"{template_number}"
            while len(template_number_str) < 4:
                template_number_str = "0" + template_number_str
            subprocess.run(["mv", f"iterations/{temp_dir}/{filename}", f"iterations/{temp_dir}/{template_number_str}.pdb"])
            print(f"###### {filename} ADDED TO {temp_dir} ({distance} A) ######")
            template_number = template_number + 1

        update_temp_dir(script_path, f"iterations/{temp_dir}")
        # Clear ouput directory
        if run_number < 2:
            clear_directory(outputdir)
        #subprocess.run(["rm", "-r", outputdir])
    return mod_count


def run_distance_finder(structure_file, p1, p2):
    """
    Runs an external script to calculate the distance between two residues
    in a protein structure.

    Args:
        structure_file (str): Path to the .pdb file.
        p1 (str): Residue index 1.
        p2 (str): Residue index 2.

    Returns:
        str or None: Distance in angstroms as a string, or None if failed.
    """
    distance = subprocess.run(
        ["python3", "distance_finder/DistanceFinder.py", structure_file, p1, p2],
        capture_output=True,
        text=True
    )
    if distance.returncode != 0:
        print("Error:", distance.stderr)
        return None
    return distance.stdout.strip()


def update_temp_dir(script_path, dir_name):
    """
    Updates the 'temp_dir' line in the shell script to point to a new directory.

    Args:
        script_path (str): Path to the shell script.
        dir_name (str): Name of the new template directory.
    """
    with open(script_path, 'r') as file:
        lines = file.readlines()
    
    with open(script_path, 'w') as file:
        for line in lines:
            if line.startswith("temp_dir="):
                file.write(f"temp_dir={dir_name}\n")
            else:
                file.write(line)


def plot_and_save_distances(distances, run_number, bin_centers):
    os.makedirs("distance_distributions", exist_ok=True)
    plot_name = f"{engine.graph_output_accuracy_bar(distances, bins=bin_centers, n=n)}"
    subprocess.run(["mv", f"{plot_name}.png", f"./distance_distributions/{plot_name}{run_number+1}.png"])
    return 0


def main():
    """
    Entry point for the ColabFold Wrapper.
    Initializes project setup and executes the template filtering loop.
    """

    # Welcome message
    print("*" * 31)
    print()
    print((" " * 7) + "ColabFold Wrapper")
    print()
    print("*" * 31)

    jobs = "jobs.json"
    mods = "mod_counts.json"

    script_path = initialize_project(jobs)

    print(">>> ATTEMPTING TO RUN COLABFOLD\n")
    n, outputdir = get_from_current_job(jobs, ["n", "outputdir"])
    n = int(n)
    mod_counts = {outputdir: {}}

    for run_number in range(5):
        """
        Start with three iterations for testing
        Once running, continue iterating until an ideal structure is output
        """
        os.chmod(script_path, 0o755)
        subprocess.run([script_path], check=True)
        iteration_mod_count = filter_output(run_number, jobs, script_path, n)
        mod_counts[outputdir][run_number] = iteration_mod_count

    # Move iterations directory into output directory to save results
    subprocess.run(["mv", "./iterations/", outputdir])
    subprocess.run(["mv", "./distance_distributions/", outputdir])

    # Append mod_counts to json logs
    append_mods_json(mods, mod_counts)


if __name__ == '__main__':
    main()
