from pathlib import Path

import remote_runner
from remote_runner import LocalWorker, Pool, Task


def main():
    remote_runner.log_to(".remote-runner.log", level="DEBUG")

    workers = [LocalWorker()]

    tasks = [
        Task.load(state)
        for state in sorted(Path("./").glob("data/output/*/*/state.dill"))
    ]

    Pool(workers).run(tasks)


if __name__ == "__main__":
    main()
