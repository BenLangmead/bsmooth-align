#!/bin/sh

[ -z "$BSMOOTH_HOME" ] && echo "BSMOOTH_HOME not set!" && exit 1
find $BSMOOTH_HOME -name '*.pm' | xargs -n 1 perl
echo "BSMOOTH_HOME was set to '$BSMOOTH_HOME'"
