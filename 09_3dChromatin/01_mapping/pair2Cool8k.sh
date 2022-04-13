#!/bin/bash

cooler cload pairix -p 6 aep.genome:8000 aep.bsorted.pairs.gz aepHic.8k.cool

#cooler zoomify --balance -p 6 aepHic.1000.cool