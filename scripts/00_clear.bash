#!/bin/bash

echo "START 00_clear"
set -e
BASE_DIR=$1

echo "Clear build and target files"
rm -rf ${BASE_DIR}/apps/*/build/*
rm -rf ${BASE_DIR}/apps/*/target/*

echo "Clear test output data"
rm -rf ${BASE_DIR}/data/testing_dynamic/*/*

echo "FINISHED 00_clear"

