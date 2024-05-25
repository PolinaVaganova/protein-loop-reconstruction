import os
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

MINUTES = 60  # seconds
DEFAULT_CHECK_CALL_TIMEOUT = 30 * MINUTES


@contextmanager
def chdir(path: Path):
    """Sets the cwd within the context

    Args:
        path (Path): The path to the cwd

    Yields:
        None
    """

    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def check_call(cmd, **kwargs) -> None:
    """
    A substitute for subprocess.check_call that suppresses output
    unless an error occurs

    :param cmd: Same as in subprocess.Popen
    :param kwargs: Same as in subprocess.Popen
    :return: None
    """

    with tempfile.TemporaryFile() as tmp_stderr:
        if "stderr" not in kwargs:
            stderr_capture = tmp_stderr
            kwargs["stderr"] = tmp_stderr
        else:
            stderr_capture = None

        if "stdout" not in kwargs:
            kwargs["stdout"] = subprocess.DEVNULL

        if "timeout" not in kwargs:
            kwargs["timeout"] = DEFAULT_CHECK_CALL_TIMEOUT

        try:
            subprocess.check_call(cmd, **kwargs)
        except (OSError, subprocess.CalledProcessError):
            if stderr_capture:
                sys.stderr.write(tmp_stderr.read().decode("utf-8"))
                raise
