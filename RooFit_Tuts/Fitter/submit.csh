#!/bin/csh

set sys = $1
set run = $2


./tmp/fitAppToy $sys 100 1 $run 
./tmp/fitAppToy $sys 100 250 $run
./tmp/fitAppToy $sys 100 500 $run
./tmp/fitAppToy $sys 100 1000 $run
./tmp/fitAppToy $sys 100 2000 $run
./tmp/fitAppToy $sys 100 4000 $run
./tmp/fitAppToy $sys 100 6000 $run

./tmp/fitAppToy $sys 120 1 $run 
./tmp/fitAppToy $sys 120 250 $run
./tmp/fitAppToy $sys 120 500 $run
./tmp/fitAppToy $sys 120 1000 $run
./tmp/fitAppToy $sys 120 2000 $run
./tmp/fitAppToy $sys 120 4000 $run
./tmp/fitAppToy $sys 120 6000 $run

./tmp/fitAppToy $sys 140 1 $run 
./tmp/fitAppToy $sys 140 250 $run
./tmp/fitAppToy $sys 140 500 $run
./tmp/fitAppToy $sys 140 1000 $run
./tmp/fitAppToy $sys 140 2000 $run
./tmp/fitAppToy $sys 140 4000 $run
./tmp/fitAppToy $sys 140 6000 $run

./tmp/fitAppToy $sys 160 1 $run 
./tmp/fitAppToy $sys 160 250 $run
./tmp/fitAppToy $sys 160 500 $run
./tmp/fitAppToy $sys 160 1000 $run
./tmp/fitAppToy $sys 160 2000 $run
./tmp/fitAppToy $sys 160 4000 $run
./tmp/fitAppToy $sys 160 6000 $run

./tmp/fitAppToy $sys 180 1 $run 
./tmp/fitAppToy $sys 180 250 $run
./tmp/fitAppToy $sys 180 500 $run
./tmp/fitAppToy $sys 180 1000 $run
./tmp/fitAppToy $sys 180 2000 $run
./tmp/fitAppToy $sys 180 4000 $run
./tmp/fitAppToy $sys 180 6000 $run
