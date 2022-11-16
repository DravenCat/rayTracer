#!/bin/sh
g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer

./RayTracer 1024 3 1 hierarchical.ppm 1

./RayTracer 1024 3 1 basic.ppm 2

./RayTracer 1024 3 1 final.ppm 3