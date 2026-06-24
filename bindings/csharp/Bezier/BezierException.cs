// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using System;

namespace Bezier;

/// <summary>Raised when the underlying C++ library reports an error (e.g. parameter out of range).</summary>
public class BezierException : Exception
{
    public BezierException(string message) : base(message) { }
}
