#!/bin/bash 
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
julia --project="${DIR}" "${DIR}/paper-2026/simulation.jl" $1 $2 $3 $4