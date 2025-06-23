import logging
from pathlib import Path
import sys


def assert_path(path: str, assert_dir: bool = True) -> Path:
    """Validates that a given path exists and is either a file or a directory.

    Converts the input string to a ``Path`` object and checks whether it exists and matches
    the expected type (file or directory), based on the ``assert_dir`` flag.

    :param str path: The file system path to validate.
    :param bool assert_dir: If ``True``, asserts that the path is a directory.
                            If ``False``, asserts that the path is a file.

    :raises TypeError: If ``path`` is not a string.
    :raises TypeError: If ``assert_dir`` is not a boolean.
    :raises FileNotFoundError: If the path does not exist.
    :raises NotADirectoryError: If ``assert_dir=True`` but the path is a file.
    :raises IsADirectoryError: If ``assert_dir=False`` but the path is a directory.

    :return: The validated path as a ``Path`` object.
    :rtype: pathlib.Path
    """
    assert isinstance(path, str), "Path input must be a string."
    assert isinstance(
        assert_dir, bool
    ), "Expected 'True' or 'False' to indicate whether a directory or file is asserted."

    # Convert path string to `Path`
    path = Path(path)

    if assert_dir:
        if path.is_dir() and path.exists():
            return path
        elif path.is_file():
            raise NotADirectoryError(f"Expected a directory, but `{path}` is a file.")
        else:
            raise FileNotFoundError(f"Directory `{path}` does not exist.")
    else:
        if path.is_file() and path.exists():
            return path
        elif path.is_dir():
            raise IsADirectoryError(f"Expected a file, but `{path}` is a directory.")
        else:
            raise FileNotFoundError(f"File `{path}` does not exist.")


def setup_logging(path_to_log: str) -> logging.Logger:
    """
    Set up and configure logging to both console and file.

    :param path_to_log: Path to the log file where logs will be written.
    :type path_to_log: str

    :returns: Configured logger instance.
    :rtype: logging.Logger
    """

    logger = logging.getLogger("logger")
    logger.setLevel(logging.DEBUG)  # Lowest log level (logs everything)

    # Custom formatter
    formatter = logging.Formatter(
        "%(asctime)s || LEVEL: %(levelname)s |> %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)  # Set handler specific log level
    console_handler.setFormatter(formatter)  # Add custom formatter to handler
    logger.addHandler(console_handler)  # Add handler to the logger

    # Add file handler
    file_handler = logging.FileHandler(path_to_log)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger
