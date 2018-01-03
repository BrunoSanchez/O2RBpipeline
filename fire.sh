#!/bin/bash
for i in {1..10}
do
    for i in {1..125}
        do
            python in_corral.py load
        done
    python in_corral.py run -s StepCrossMatch
    python in_corral.py run -s StepSCrossMatch
    python in_corral.py run -s StepSCorrCrossMatch
    python in_corral.py run -s StepCrossMatchOIS
    python in_corral.py run -s StepCrossMatchHOT
done

