import argparse
from pele.concurrent import ConnectWorker
from start_server import create_system, get_server_uri, get_database_params_worker

def main():
    parser = argparse.ArgumentParser(description="connect worker queue")
    parser.add_argument("p", type=int, help="p-spin")
    parser.add_argument("nspins", type=int, help="number of spins")
    parser.add_argument("--strategy", type=str, help="strategy to adopt: random (default), "
                                                     "untrap, combine, gmin", default="random")
    args = parser.parse_args()

    nspins = args.nspins
    p = args.p

    interactions = get_database_params_worker(nspins, p)
    system = create_system(nspins, p, interactions)

    uri = get_server_uri(nspins, p)
    worker = ConnectWorker(uri, system=system, strategy=args.strategy)
    worker.run()


if __name__ == "__main__":
    main()
