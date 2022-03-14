#!/bin/bash

(cd build; make && cd .. && rm -f TestFiles/SW*CH7*.cdf && ./build/bin/chaos A20220309 shc TestFiles TestFiles)

