import pathlib

def get_include() -> pathlib.Path:
    import pele
    include_path = pathlib.Path(pele.__file__).parent / "include"
    current_path = pathlib.Path(".").absolute()
    if include_path.is_relative_to(current_path):
        return include_path.relative_to(current_path)
    else:
        return include_path
