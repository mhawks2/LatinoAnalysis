#!/bin/bash

github-addext() {

    git clone --quiet git@github.com:$1 $2
    if [ ! -z $3 ]; then
    (
        cd $2
        git checkout --quiet $3
    )
    fi

}
