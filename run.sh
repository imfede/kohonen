#!/bin/bash

rm -rf datas/*
rm *.class

make -k

java Kohonen

gnuplot graph.plt
