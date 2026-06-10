import subprocess
import os

r_script = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\PC_30genes_analysis.R"
r_exe = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\AI 代码编写\R-4.5.2\bin\x64\Rscript.exe"

# Run the R script
result = subprocess.run([r_exe, r_script], capture_output=True, text=True)
print("STDOUT:")
print(result.stdout)
print("\nSTDERR:")
print(result.stderr)
print("\nReturn code:", result.returncode)
