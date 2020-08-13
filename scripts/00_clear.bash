#!/bin/bash

echo "Clear build and dist files"
rm -r ../apps/*/build/*
rm -r ../apps/*/dist/*

echo "Clear test output data"
rm -r ../data/testing_dynamic/*/*

echo "ls build and dist files"
ls ../apps/*/build
ls ../apps/*/dist

echo "done"
