#!/bin/bash
git config --global user.name "Nils Huesken"
git config --global user.email nils.hueskelgooglemail.com
git add .
git commit -m "$1"
git push origin userConstraints
