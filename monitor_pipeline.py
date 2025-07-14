import subprocess
import sys
from datetime import datetime


def main():
    log_path = 'pipeline_output.log'
    with open(log_path, 'w') as f:
        process = subprocess.Popen(
            [sys.executable, 'interactive_trophoblast_pipeline.py'] + sys.argv[1:],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        for line in process.stdout:
            timestamped = f"{datetime.now().isoformat()} | {line}"
            print(timestamped, end='')
            f.write(timestamped)
        process.wait()

if __name__ == '__main__':
    main()
