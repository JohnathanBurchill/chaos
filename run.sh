#!/bin/bash

(cd build; make && cd .. && ./build/bin/chaos A20220309 shc TestFiles TestFiles)

