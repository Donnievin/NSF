import subprocess
import os

# Define the list of scripts to be executed in order
scripts = [
    "Pipeline/0_PRSM.py",
    "Pipeline/1_CLEAN.py",
    "Pipeline/2_MASK.py",
    "Pipeline/3_ESM3.py",
    "analysis/1_Dim_Red.py",
    "analysis/2_Properties.py",
    "analysis/3_RMSD.py",
    "analysis/9_Fast_Graphs.py"
]

# Function to execute a script
def execute_script(script):
    try:
        print(f"Executing {script}...")
        subprocess.run(["python3", script], check=True)
        print(f"{script} executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing {script}: {e}")
        exit(1)

# Main function to run the scripts in order
def main():
    # Creates an environment from the yaml file
    try:
        os.system('conda env create -f environment.yml')
    except Exception as e:
        print("Could not find the environment file")

    # Activates the environment 
    os.system('conda init')
    os.system('conda activate nsf_env')

    # Runs the codes in order as shown in the list above

    for script in scripts:
        if os.path.exists(script):
            execute_script(script)
        else:
            print(f"Script {script} does not exist.")

if __name__ == "__main__":
    main()
