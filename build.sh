#!/usr/bin/env bash
mkdir -p ./build
rm -rf ./build/*
g++ -O -Wall -std=c++20 -shared -fpic -o ./build/info.so ./cpp/*.cpp