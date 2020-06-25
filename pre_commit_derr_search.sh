#!/usr/bin/env bash
FORBIDDEN='derr.*<<'
git diff | \
	grep $FORBIDDEN | \
	echo "COMMIT REJECTED Found '$FORBIDDEN' references. Please remove them before commiting" && exit 1