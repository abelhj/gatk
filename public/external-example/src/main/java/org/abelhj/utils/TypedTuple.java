package org.abelhj.utils;

import java.io.Serializable;

public class TypedTuple<L, R> implements Serializable {
    private static final long serialVersionUID = 1L;

    protected L left;
    protected R right;

    protected TypedTuple() {
	// Default constructor for serialization
    }

    public TypedTuple(L inLeft, R inRight) {
	left = inLeft;
	right = inRight;
    }

    public L getLeft() {
	return left;
    }

    public R getRight() {
	return right;
    }
}