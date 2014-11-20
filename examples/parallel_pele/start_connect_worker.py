from pele.concurrent import ConnectWorker

from start_server import create_system, get_server_uri


def main():
    system = create_system()

    uri = get_server_uri()
    worker = ConnectWorker(uri, system=system, strategy="random")
    worker.run()


if __name__ == "__main__":
    main()
