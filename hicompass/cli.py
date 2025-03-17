#!/usr/bin/env python3
import argparse
import sys
from hicompass.commands import predict

def main():
    """
    Hi-Compass main command
    """
    parser = argparse.ArgumentParser(
        description='Hi-Compass: A Hi-C predict deep learning model',
        usage='hicompass <command> [<args>]'
    )
    
    # 添加子命令解析器
    subparsers = parser.add_subparsers(dest='command', help='commands that you can use')
    
    # 添加predict子命令
    predict_parser = subparsers.add_parser('predict', help='predict a Hi-C cool result from input data')
    predict.configure_parser(predict_parser)
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 如果没有指定子命令，显示帮助信息
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    # 根据子命令执行相应的功能
    if args.command == 'predict':
        predict.run(args)

if __name__ == "__main__":
    main()