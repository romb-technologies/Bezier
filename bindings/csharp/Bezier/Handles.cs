// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using System;
using System.Runtime.InteropServices;

namespace Bezier;

internal sealed class CurveHandle : SafeHandle
{
    public CurveHandle() : base(IntPtr.Zero, ownsHandle: true) { }

    public CurveHandle(IntPtr h) : base(IntPtr.Zero, ownsHandle: true) => SetHandle(h);

    public override bool IsInvalid => handle == IntPtr.Zero;

    protected override bool ReleaseHandle()
    {
        Native.bz_curve_destroy(handle);
        return true;
    }
}

internal sealed class PolyCurveHandle : SafeHandle
{
    public PolyCurveHandle() : base(IntPtr.Zero, ownsHandle: true) { }

    public PolyCurveHandle(IntPtr h) : base(IntPtr.Zero, ownsHandle: true) => SetHandle(h);

    public override bool IsInvalid => handle == IntPtr.Zero;

    protected override bool ReleaseHandle()
    {
        Native.bz_polycurve_destroy(handle);
        return true;
    }
}
