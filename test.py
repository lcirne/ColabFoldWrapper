import argparse

if __name__=='__main__':
    parser  = argparse.ArgumentParser()
    parser.add_argument("--test", "-t",
            default="1")
    args = parser.parse_args()
    print(args.test)
