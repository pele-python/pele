from pele.concurrent import BasinhoppingWorker

from start_server import create_system, get_server_uri


def main():
    system = create_system()

    uri = get_server_uri()
    worker = BasinhoppingWorker(uri, system=system)
    worker.run(1000)


if __name__ == "__main__":
    main()
