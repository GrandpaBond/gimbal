/*
This code will form the basis for the extension "heading.pxt".
The top-level functionality will then be moved into a separate "test.ts" module.

*/
// =============  GLOBALS ===============

const Window = 7
const SampleGap = 20
const EnoughScanTime = 1500 // minimum acceptable scan-time
const EnoughSamples = 70 // fewest acceptable scan samples
const TooManySamples = 500 // don't be too greedy with memory!
const MarginalField = 10 // minimum acceptable field-strength for magnetometer readings
const TinyField = 2 // minimal field magnitude, considered to be a zero-crossing

enum Tasks {
    PerformScan,
    SetNorth,
    TakeBearing,
}
let nextTask: Tasks = Tasks.PerformScan

let simulating = isSimulating() // true when debugging

let scan: Scan  // array of scanned magnetometer samples
let testReadings: Reading[] = [] // array of test readings
let testIndex: number // next testReading to use when simulating

// Sensor Measurements
let magnetXYZ: Vector // current magnetic field
let gravityXYZ: Vector // current accelerometer pose
let startXYZ: Reading // reading of starting field and pose of the buggy (deemed north and upright)

let northXYZ: Vector // starting magnetic field of the buggy (while pointing "North")
let downXYZ: Vector // buggy's Down axis measured while upright and stationary(fixed, dependent on mounting)

// calibrated correction adjustments for accelerometer readings (adopting explicit calibration values
// previously measured externally by reading the extreme static values in each dimension)
let poseOffset: Vector // central offsets from origin in each dimension
let poseScaleY: number // multiplier to match Y readings with X
let poseScaleZ: number // multiplier to match Z readings with X



// re-orientation rotations
let rotateXYZtoRFD: Quaternion // sensor [XYZ] to buggy's [Right,Front,Down] frame 
let rotateRFDtoENG: Quaternion // buggy [Right,Front,Down] to world [East,North,Gravity] frame 
let rotateXYZtoENG: Quaternion // sensor [XYZ] directly to world [East,North,Gravity] frame

let magnetENG: Vector
let gravityENG: Vector
let heading: number

// ================ CLASSES ================

/* 3-D vector, with methods for normalisation, dot-product and cross-product. 
*/
class Vector {
    x: number
    y: number
    z: number

    constructor(dx: number, dy: number, dz: number) {
        this.x = dx
        this.y = dy
        this.z = dz
    }

    normalised(): Vector {
        let r = this.getMagnitude()
        if (r == 0) {
            return new Vector(0, 0, 0)
        } else {
            return new Vector(this.x / r, this.y / r, this.z / r)
        }
    }

    dottedWith(v: Vector): number {
        return (this.x * v.x + this.y * v.y + this.z * v.z)
    }

    crossedWith(v: Vector): Vector {
        let x = this.y * v.z - v.y * this.z
        let y = this.z * v.x - v.z * this.x
        let z = this.x * v.y - v.x * this.y
        return new Vector(x, y, z)
    }

    // we are sometimes more interested in the square of the magnitude 
    getLengthSquared(): number {
        return ((this.x * this.x) + (this.y * this.y) + (this.z * this.z))
    }

    getMagnitude(): number {
        return Math.sqrt(this.getLengthSquared())
    }
}

/**
* A Quaternion is a tool for manipulating rotations.
* Initial construction is from an rotation about a given axis.
* Tools are provided to make it represent an alignment between two vectors,
* and to apply it to rotate a vector.
*/
class Quaternion {
    // the real part
    w: number
    // the three imaginary parts
    i: number
    j: number
    k: number
    // squares of components (precomputed for efficiency)
    ww: number
    ii: number
    jj: number
    kk: number
    // doubled products of components (precomputed for efficiency)
    wi2: number
    wj2: number
    wk2: number
    ij2: number
    jk2: number
    ki2: number

    // given a rotation-angle and an axis-direction, build a unit quaternion
    constructor(angle: number, axis: Vector) {
        let unitV = axis.normalised()
        this.w = Math.cos(angle / 2)
        let sinHalfAngle = Math.sin(angle / 2)
        this.i = unitV.x * sinHalfAngle
        this.j = unitV.y * sinHalfAngle
        this.k = unitV.z * sinHalfAngle
        this.precompute()
    }

    // compute the Quaternion needed to align unit vector (a) onto unit vector (b)
    // by rotating about an axis normal to their common plane
    // The axis is just their cross-product, and the angle is deived from their dot-product.
    // We need to deal with two special cases:
    //  - if (a) parallel to (b) (or very close), use the unit Quaternion (1,0,0,0)
    //  - if (a) opposite to (b) (or very close), we need a rotation by 180 degrees around any orthogonal axis

    toAlignVectors(a: Vector, b: Vector) {
        let cross: Vector
        let dot = a.dottedWith(b)
        cross = a.crossedWith(b) // degenerates if vectors align 
        if (dot > 0.999999) {
            // parallel vectors: use identity quaternion 
            this.w = 1
            this.i = 0
            this.j = 0
            this.k = 0
        } else if (dot < -0.999999) {
            // anti-parallel vectors need 180-degree rotation about any orthogonal axis.
            // first try using the normal to the plane through (a) and the x-axis
            cross = a.crossedWith(new Vector(1, 0, 0))
            // if by chance (a) WAS the x-axis, use the y-axis as our rotation axis
            if (cross.getLengthSquared() < 0.00001) {
                cross = new Vector(0, 1, 0)
            }
            this.i = cross.x
            this.j = cross.y
            this.k = cross.z
            // the half-angle needed will be 90 degrees
            this.w = Math.PI
        } else {
            // general case
            this.i = cross.x
            this.j = cross.y
            this.k = cross.z
            this.w = dot
        }
        this.normalise()
        this.precompute()

        datalogger.log(datalogger.createCV("w", this.w),
            datalogger.createCV("i", this.i),
            datalogger.createCV("j", this.j),
            datalogger.createCV("k", this.k))
    }



    // use this Quaternion to generate a rotated Vector
    appliedToVector(v: Vector): Vector {
        let result = new Vector(0, 0, 0)
        result.x
            = v.x * (this.ww + this.ii - this.jj - this.kk)
            + v.y * (this.ij2 - this.wk2)
            + v.z * (this.ki2 + this.wj2)

        result.y
            = v.y * (this.ww + this.jj - this.kk - this.ii)
            + v.z * (this.jk2 - this.wi2)
            + v.x * (this.ij2 + this.wk2)

        result.z
            = v.z * (this.ww + this.kk - this.ii - this.jj)
            + v.x * (this.ki2 - this.wj2)
            + v.y * (this.wi2 + this.jk2)

        return result
    }

    // for a unit Quaternion, the squares of all its components add up to 1.
    normalise() {
        this.ww = this.w * this.w
        this.ii = this.i * this.i
        this.jj = this.j * this.j
        this.kk = this.k * this.k
        let r = Math.sqrt(this.ww + this.ii + this.jj + this.kk)
        this.w /= r
        this.i /= r
        this.j /= r
        this.k /= r
    }


    // precompute squares and products (some doubled)...
    precompute() {
        this.ww = this.w * this.w
        this.ii = this.i * this.i
        this.jj = this.j * this.j
        this.kk = this.k * this.k
        this.wi2 = this.w * this.i * 2
        this.wj2 = this.w * this.j * 2
        this.wk2 = this.w * this.k * 2
        this.ij2 = this.i * this.j * 2
        this.ki2 = this.i * this.k * 2
        this.jk2 = this.j * this.k * 2
    }
}

// a Reading is a compound object containing synchronous 3-D readings from both the magnetometer and accelerometer
class Reading {
    field: Vector // average magnetometer reading
    pose: Vector // average accelerometer reading

    constructor(fieldX: number, fieldY: number, fieldZ: number,
        poseX: number, poseY: number, poseZ: number) {
        this.field = new Vector(fieldX, fieldY, fieldZ)
        this.pose = new Vector(poseX, poseY, poseZ)
    }
}

/* a Sample is a time-stamped 3-D reading from the magnetometer (one element of the scan array)
class Sample {
    time: number
    field: Vector

    constructor(t: number, fieldX: number, fieldY: number, fieldZ: number) {
        this.time = t
        this.field = new Vector(fieldX, fieldY, fieldZ)
    }
}*/


/* A Smoother object computes moving averages from a sequence of time-stamped vectors of values.
    It is used to smooth out jittery sensors such as the magnetometer or accelerometer.
    Timing irregularites due to scheduler interrupts demand this somewhat complex maths.
    The number of readings to be averaged (window) and expected gap between readings (samplingGap)
    together set the overall latency associated with the exponential averaging process
    and govern the blending of new and old readings.
    */

class Smoother {
    dims: number; // dimensionality
    averages: number[] = []; // the rolling averages
    window: number; // number of samples needed to form a good average
    samplingGap: number; // time gap between expected readings
    latency: number // resulting time taken to collect a good moving average from scratch
    lastTime: number; // timestamp of latest readings
    lastInputs: number[] = []; // copy of latest set of readings

    constructor(startTime: number, window: number, samplingGap: number, initialValues: number[]) {
        this.reset(startTime, window, samplingGap, initialValues)
    }

    // (re)initialise this Smoother
    reset(startTime: number, window: number, samplingGap: number, initialValues: number[]) {
        this.lastTime = startTime
        this.window = window
        this.samplingGap = samplingGap
        this.latency = window * samplingGap
        this.dims = initialValues.length
        for (let dim = 0; dim < this.dims; dim++) {
            this.averages[dim] = initialValues[dim]
            this.lastInputs[dim] = initialValues[dim]
        }
    }

    update(timeStamp: number, values: number[]): number[] {
        // work out appropriate blend, based on time-step (guarding against zero!)
        let timeFraction = (timeStamp - this.lastTime + 1) / this.latency
        let keepOld = Math.exp(-timeFraction)
        let inherited = (1 - keepOld) / timeFraction
        // amplify the most recent sample's contribution to the inherited average
        let boostLast = (inherited - keepOld)
        let addNew = (1 - inherited)
        // (blending proportions keepOld + boostLast + addNew will always add up to 100%)
        // apply blending to all elements of old and new data arrays
        let result: number[] = []
        for (let i = 0; i < this.dims; i++) {
            result.push((keepOld * this.averages[i])
                + (boostLast * this.lastInputs[i])
                + (addNew * values[i]))
        }
        // update history for next time around
        this.averages = result
        this.lastTime = timeStamp
        this.lastInputs = values

        return result
    }
}

/** A Scan is a dataset of sequential magnetometer readings gathered while the buggy is spinning on the spot.
 * Methods are provided to acquire, scope and analyse this sequence to derive the correction parameters
 * for the magnetometer (used for future readings). 
 * Analysis of the dataset also reveals how long each rotation took, and the orientation of the spin-axis 
 * (measured in the sensor's XYZ frame).
 * 
*/
class Scan {
    samples: Vector[] // sequence of magnetometer & accelerometer readings
    times: number[] // matching sequence of time-stamps for fields[]
    period: number // derived spin-rotation period in ms
    downXYZ: Vector // spin-axis (giving the buggy's "Down" axis in sensor coordinates)
    range: Vector   // field amplitudes in each dimension
    strength: number // the average magnetic field-strength detected on a scan 

    // calibrated correction adjustments for magnetometer readings
    fieldOffset: Vector  // central offsets from origin in each dimension
    fieldScaleY: number // multiplier to match Y readings with X
    fieldScaleZ: number // multiplier to match Z readings with X

    fieldSmoother: Smoother // uses a Smoother to maintain a rolling average
    constructor() {
        this.samples = []
        this.times = []
    }

    // SCAN METHODS


    // Perform a scan for specified time
    acquire(ms: number, dumpIt: boolean) {
        let timeWas: number
        let timeNow: number
        this.samples = [] // start with empty array
        this.times = []

        // get initial reading
        let timeStamp = input.runningTime()
        let field: number[] = [
            input.magneticForce(Dimension.X),
            input.magneticForce(Dimension.Y),
            input.magneticForce(Dimension.Z)]

        this.fieldSmoother = new Smoother(timeStamp, Window, SampleGap, field)
        let smooth: number[]

        // after an initial settling period, continue cranking out updated moving averages... 
        let startTime = timeStamp + (Window * SampleGap)
        let stopTime = timeStamp + ms

        // ...until we run out of time (or space!)
        while ((timeStamp < stopTime)
            && (this.samples.length < TooManySamples)) {
            // After processing, sleep until it's time for next sample.
            // NOTE: here is where various system subprograms will get scheduled.
            // If they need more time than we've offered, our next sample will get delayed!
            // (This seems to incur extra delays of ~44 ms every 100ms, plus ~26ms every 400ms)

            timeWas = timeStamp // remember time of latest sample
            timeNow = input.runningTime()
            basic.pause((timeWas + SampleGap) - timeNow) // pause for remainder of SampleGap (if any!)
            timeStamp = input.runningTime() // take a fresh set of readings

            field = [
                input.magneticForce(Dimension.X),
                input.magneticForce(Dimension.Y),
                input.magneticForce(Dimension.Z)]
            smooth = this.fieldSmoother.update(timeNow, field)

            // only start recording once the moving average has stabilised
            if (timeStamp > startTime) {
                // store the averaged field values (as a deep copy!)
                this.samples.push(new Vector(smooth[0], smooth[1], smooth[2]))
                this.times.push(timeNow)  // timestamp it  
            }
        }

        // dump this scan to the datalogger
        if (dumpIt) {
            for (let i = 0; i < this.samples.length; i++) {
                datalogger.log(
                    datalogger.createCV("data", "raw scan"),
                    datalogger.createCV("fx", this.samples[i].x),
                    datalogger.createCV("fy", this.samples[i].y),
                    datalogger.createCV("fz", this.samples[i].z))
            }
        }
    }


    // Each dimension should track a sinusoidal wave of values (generally not centred on zero).
    // This method finds the value ranges for each axis (usually NOT the full field-strength in any dimension)
    // It also sets the global offsets needed to correctly re-centre biased future readings
    scope() {
        let xlo = 9999999
        let ylo = 9999999
        let zlo = 9999999
        let xhi = -9999999
        let yhi = -9999999
        let zhi = -9999999
        for (let i = 0; i < this.samples.length; i++) {
            xhi = Math.max(xhi, this.samples[i].x)
            yhi = Math.max(yhi, this.samples[i].y)
            zhi = Math.max(zhi, this.samples[i].z)
            xlo = Math.min(xlo, this.samples[i].x)
            ylo = Math.min(ylo, this.samples[i].y)
            zlo = Math.min(zlo, this.samples[i].z)
        }

        // derive RMS field-strength from the ranges detected in each axis
        let rangeX = (xhi - xlo) / 2
        let rangeY = (yhi - ylo) / 2
        let rangeZ = (zhi - zlo) / 2
        this.range = new Vector(rangeX, rangeY, rangeZ)
        this.strength = Math.sqrt((rangeX * rangeX) + (rangeY * rangeY) + (rangeZ * rangeZ))

        // offsets from the origin (due to "hard-iron" distortions) lie mid-way between extremes
        let offX = (xhi + xlo) / 2
        let offY = (yhi + ylo) / 2
        let offZ = (zhi + zlo) / 2
        this.fieldOffset = new Vector(offX, offY, offZ)
    }

    recentre() {
        // re-centre all the scan samples, eliminating "hard-iron" environmental magnetic effects.
        for (let i = 0; i < this.samples.length; i++) {
            this.samples[i].x -= this.fieldOffset.x
            this.samples[i].y -= this.fieldOffset.y
            this.samples[i].z -= this.fieldOffset.z
        }
    }

    // Method to analyse the scan-readings and derive the magnetometer scaling factors
    // and the scan spin-axis (measured in the XYZ sensor frame).
    analyse() {
        /* given the set of six [X,Y,Z] measurements:
                [M, N, -] when crossing the XY plane
                [-, P, Q] when crossing the YZ plane
                [R, -, S] when crossing the ZX plane
    
        ...and knowing that: 
                X**2 + (yScale * Y)**2 + (zScale * Z)**2 = B**2 (the square of the field strength)
        
        ...we can (after some maths!) derive the calibration factors (relative to x):
                yScale = sqrt((MMQQ - MMSS - QQRR) / (SSNN - SSPP - NNQQ))
                zScale = sqrt((PPRR - PPMM - RRNN) / (SSNN - SSPP - NNQQ))
        */

        // we'll mostly be using the squares of the zero-crossing components
        let MM = 0
        let NN = 0
        let PP = 0
        let QQ = 0
        let RR = 0
        let SS = 0
        // preserve history
        let xWas: number
        let yWas: number
        let zWas: number

        // First, collect the plane-crossings in each direction.
        // Simultaneously, collect half-periods of rotation, which we will average.

        // counts of zero-crossings detected in this scan
        let nCrossXY = 0
        let nCrossYZ = 0
        let nCrossZX = 0
        // time-stamps of first crossings (not yet found)
        let xStart = -1
        let yStart = -1
        let zStart = -1
        // timestamps of last crossings
        let xFinish: number
        let yFinish: number
        let zFinish: number

        // flags to inhibit clocking multiple jittery crossings 
        let needXY = true
        let needYZ = true
        let needZX = true

        let x = this.samples[0].x
        let y = this.samples[0].y
        let z = this.samples[0].z
        
        for (let i = 0; i < this.samples.length; i++) {
            xWas = x
            yWas = y
            zWas = z
            x = this.samples[i].x
            y = this.samples[i].y
            z = this.samples[i].z

            // avoid any exact zeroes (they only complicate comparisons!)
            if (x == 0) x = xWas
            if (y == 0) y = yWas
            if (z == 0) z = zWas

            // Look for the first transition of each half-cycle (i.e. where the sign flips)
            // (jitter or near-axis alignment may cause repeated fluctuations, which we ignore)

            if ((z * zWas < 0) && needXY) { // sign of z value flips when crossing the XY plane
                MM += x ** 2
                NN += y ** 2
                nCrossXY++
                zFinish = this.times[i]
                if (zStart < 0) zStart = zFinish // start the clock...
                needXY = false
                // got this plane-crossing, so now only allow other planes to be detected
                needYZ = true
                needZX = true
            }
            if ((x * xWas < 0) && needYZ) { // sign of x value flips when crossing the YZ plane
                PP += y ** 2
                QQ += z ** 2
                nCrossYZ++
                xFinish = this.times[i]
                if (xStart < 0) xStart = xFinish
                needYZ = false
                needXY = true
                needZX = true
            }
            if ((y * yWas < 0) && needZX) { // sign of y value flips when crossing the ZX plane
                RR += x ** 2
                SS += z ** 2
                nCrossZX++
                yFinish = this.times[i]
                if (yStart < 0) yStart = yFinish
                needZX = false
                needXY = true
                needYZ = true
            }
        }
        // average the squared crossing points
        MM /= nCrossXY
        NN /= nCrossXY
        PP /= nCrossYZ
        QQ /= nCrossYZ
        RR /= nCrossZX
        SS /= nCrossZX

        // derive the average "flip" times (each making half a rotation)
        let xFlip = (xFinish - xStart) / (nCrossYZ - 1)
        let yFlip = (yFinish - yStart) / (nCrossZX - 1)
        let zFlip = (zFinish - zStart) / (nCrossXY - 1)

        // average the three half-periods, then double them to get our best measure for full period
        this.period = (xFlip + yFlip + zFlip) / 1.5

        // construct the relative scaling factors
        let bottom = (NN * SS) - (SS * PP) - (NN * QQ)
        this.fieldScaleY = Math.sqrt((MM * QQ) - (QQ * RR) - (SS * MM) / bottom)
        this.fieldScaleZ = Math.sqrt((PP * RR) - (PP * MM) - (NN * RR) / bottom)

        /* retrospectively rebalance the Y and Z components of the plane-crossing vectors
                [M, N, -] when crossing the XY plane
                [-, P, Q] when crossing the YZ plane
                [R, -, S] when crossing the ZX plane
        */
        let M = Math.sqrt(MM)
        let N = Math.sqrt(NN) * this.fieldScaleY
        let P = Math.sqrt(PP) * this.fieldScaleY
        let Q = Math.sqrt(QQ) * this.fieldScaleZ
        let R = Math.sqrt(RR)
        let S = Math.sqrt(MM) * this.fieldScaleZ

        // Since the three crossing-points form a co-planar triangle lying in the Spin-Circle plane, we can take the 
        // cross-product of any two edges to derive dynamically the orthogonal rotation-axis (the buggy's "Down" axis).
        // (We'll later compare this with the static reading taken when setNorth() is invoked.)
        let I = (Q * N) - (N * S) + (S * P)
        let J = (R * Q) - (Q * M) + (M * S)
        let K = (N * R) - (R * P) + (P * M)

        this.downXYZ = new Vector(I, J, K)
        this.downXYZ = this.downXYZ.normalised()

        let check = 0 // just a debug point...
    }

    /* adopt a previously-recorded dataset
    use(samples: Vector[], times: number[]) {
        this.samples = samples
        this.times = times
    }*/

    // dump the correction parameters and spin-axis
    dumpAnalysis() {
        datalogger.log(
            datalogger.createCV("yScale", this.fieldScaleY),
            datalogger.createCV("zScale", this.fieldScaleZ),
            datalogger.createCV("downX", this.downXYZ.x),
            datalogger.createCV("downY", this.downXYZ.y),
            datalogger.createCV("downZ", this.downXYZ.z))
    }
    
}


// ============== INPUT HANDLERS ===============
input.onButtonPressed(Button.A, function() {
    doNextTask()
})
input.onButtonPressed(Button.B, function () {
    dumpTestData()
})

input.onButtonPressed(Button.AB, function () {
    datalogger.deleteLog()
    basic.showIcon(IconNames.No)
    pause(2000)
    basic.clearScreen()
    nextTask = Tasks.PerformScan
    characteriseAccelerometer() // adopt calibration data for well-known (to me!) microbits
})




/**
     * Although fairly close, the magnetometer sensitivity in each axis direction varies by a few
     * percent. By extracting plane-crossings from the scan-data this function calculates from first
     * principles the global calibration factors: yScale and zScale.
     * These are then used to correct the plane-crossings before using them to derive the spin-axis.
     * As a by-product, the sample timestamps allow the average spin-rotation period to be measured.
     *
     * NOTE: There is no guarantee that the spin-axis is truly "vertical": the buggy may be operating
     * on a tilted surface. Its "Down" axis would not then coincide with the world-frame "Gravity" axis.
     * To establish this relationship, we will need (later) to call SetNorth() with the buggy at rest.
    */

// ============== FUNCTIONS ===============



function doNextTask() {
    let bearing: number
    let result: number
    switch (nextTask) {
        case Tasks.PerformScan:
            basic.showString("S") // scan
            pause(1000)
            basic.clearScreen()
            if (isSimulating) {
                result = simulateScan("T07260757_dash70")
            } else {
                scan.acquire(6000, true)
            }
            
            scan.scope() // find extremes of rotational variation
            // TODO. check here that scan.strength is sufficient

            scan.recentre() // correct for "hard-iron" bias

            scan.analyse()  // derive rotation-period and rotation-axis
            result = 0
            if (result != 0) {
                basic.showNumber(result)
            } else {
                scan.samples = [] // release memory used for scan data...
                scan.times = [] // .. and their timestamps
                basic.showIcon(IconNames.Yes)
                pause(1000)
                nextTask = Tasks.SetNorth
            }
            break

        case Tasks.SetNorth:
            basic.showString("N")
            pause(500) // ensure accelerometer is at rest
            setNorth()  // take a fix on "North" and the "Down" orientation
            pause(1000)
            basic.clearScreen()
            nextTask = Tasks.TakeBearing
            break

        case Tasks.TakeBearing:
            bearing = getHeading()
            basic.showNumber(bearing)
            pause(1000)
            basic.clearScreen()
            break

    }
}

/***
 * function correctedField(): Vector {
    let reading = new Vector(0, 0, 0)
    if (simulating) {
        reading.x = 8.16
        reading.y = 7.91
        reading.z = 32.72
    } else {
        reading.x = (input.magneticForce(0) - fxOff)
        reading.y = (input.magneticForce(1) - fyOff) * fyScale
        reading.z = (input.magneticForce(2) - fzOff) * fzScale
    }
    return reading
}

function correctedGravity(): Vector {
    let reading = new Vector(0, 0, 0)
    if (simulating) {
        reading.x = -23.53
        reading.y = 30.43
        reading.z = -762.48
    } else {
        reading.x = (input.acceleration(0) - poseOffset.x)
        reading.y = (input.acceleration(1) - poseOffset.y) * gyScale
        reading.z = (input.acceleration(2) - poseOffset.z) * gzScale
    }
    return reading
}
***/

// either we're simulating, or we're shut in a magnetic shielding box!
function isSimulating(): boolean {
    let x = input.magneticForce(0)
    let y = input.magneticForce(1)
    let z = input.magneticForce(2)
    return ((x == 0) && (y == 0) && (z == 0))
}

/* eventual user interfaces

function scanClockwise(ms: number): number {

    let nSamples = scan.samples.length

    // Now analyse the scan-data to decide how best to use the magnetometer readings.
    // we'll typically need about a couple of second's worth of scanned readings...
    let scanDuration = scan.times[scan.samples.length - 1] = scan.times[0]
    if ((this.samples.length < EnoughSamples) || (scanDuration < EnoughScanTime)) {
        return -1 // "NOT ENOUGH SCAN DATA"
    }

    let strength = scan.scope()

    // Complain if the scan didn't properly detect the Earth's magnetic field,
    // (perhaps due to magnetic shielding?)
    if (strength < MarginalField) {
        return -2 // "FIELD STRENGTH TOO WEAK"
    }
}


    // assess the scan-data to detect unequal axis sensitivity 
    // (also derives the scanPeriod, and the downXYZ spin-axis)
    // analyseScan()

    /* correct all the scan-data (for unequal axis sensitivity) by rescaling y & z values
    for (let i = 0; i < this.samples.length; i++) {
        scan[i].field.y *= yScale
        scan[i].field.z *= zScale
    }
*/


function getHeading() {
    let reading: Reading = takeReading()
    magnetXYZ = reading.field
    gravityXYZ = reading.pose
    datalogger.log(
        datalogger.createCV("data", "XYZ vals"),
        datalogger.createCV("fx", magnetXYZ.x),
        datalogger.createCV("fy", magnetXYZ.y),
        datalogger.createCV("fz", magnetXYZ.z),
        datalogger.createCV("gx", gravityXYZ.x),
        datalogger.createCV("gy", gravityXYZ.y),
        datalogger.createCV("gz", gravityXYZ.z))
    //let dot = field.dottedWith(gravity)
    //let cross = field.crossedWith(gravity)
    magnetENG = rotateXYZtoENG.appliedToVector(magnetXYZ)
    gravityENG = rotateXYZtoENG.appliedToVector(gravityXYZ)

    datalogger.log(
        datalogger.createCV("data", "ENG vals"),
        datalogger.createCV("fx", magnetENG.x),
        datalogger.createCV("fy", magnetENG.y),
        datalogger.createCV("fz", magnetENG.z),
        datalogger.createCV("gx", gravityENG.x),
        datalogger.createCV("gy", gravityENG.y),
        datalogger.createCV("gz", gravityENG.z))

    heading = (2 * Math.PI + Math.atan2(magnetENG.y, magnetENG.x)) % (2 * Math.PI)
    heading = heading * 180 / Math.PI
    datalogger.log(
        datalogger.createCV("heading", heading))
    return heading
}


// dump the test readings from this session to the datalogger
function dumpTestData() {
    for (let i = 0; i < testReadings.length; i++) {
        datalogger.log(
            datalogger.createCV("data", "test readings"),
            datalogger.createCV("fx", testReadings[i].field.x),
            datalogger.createCV("fy", testReadings[i].field.y),
            datalogger.createCV("fz", testReadings[i].field.z),
            datalogger.createCV("gx", testReadings[i].pose.x),
            datalogger.createCV("gy", testReadings[i].pose.y),
            datalogger.createCV("gz", testReadings[i].pose.z))

    }
}

// take (stable!) sensor readings for buggy "Down" axis pose and "North" magnetic field
// (measured in the sensor's XYZ frame)
function setNorth() {
    let reading: Reading
    if (simulating) {
        reading = testReadings[testIndex]
        testIndex++
    } else {
        reading = takeReading()
    }
    northXYZ = new Vector(reading.field.x, reading.field.y, reading.field.z)
    downXYZ = new Vector(reading.pose.x, reading.pose.y, reading.pose.z)

    datalogger.log(
        datalogger.createCV("data", "N & DOWN"),
        datalogger.createCV("fx", northXYZ.x),
        datalogger.createCV("fy", northXYZ.y),
        datalogger.createCV("fz", northXYZ.z),
        datalogger.createCV("gx", downXYZ.x),
        datalogger.createCV("gy", downXYZ.y),
        datalogger.createCV("gz", downXYZ.z))

 // compute rotation required to convert XYZ readings into the East-North-Gravity world-frame
    let vertical = new Vector(0, 0, -1023)
    rotateXYZtoENG = new Quaternion(0,vertical)
    rotateXYZtoENG.toAlignVectors(downXYZ, vertical)
}



// take a single test reading in the XYZ sensor-frame
function takeReading(): Reading {
    let reading: Reading
    // field accumulator
    let fieldX: number
    let fieldY: number
    let fieldZ: number
    // pose accumulator
    let poseX: number
    let poseY: number
    let poseZ: number
    if (simulating) {
        reading = testReadings[testIndex]
        testIndex++
    } else {
        for (let i = 0; i < Window; i++) {
            fieldX += input.magneticForce(Dimension.X)
            fieldY += input.magneticForce(Dimension.Y)
            fieldZ += input.magneticForce(Dimension.Z)
            poseX += input.acceleration(Dimension.X)
            poseY += input.acceleration(Dimension.Y)
            poseZ += input.acceleration(Dimension.Z)
        }
        fieldX /= Window
        fieldY /= Window
        fieldZ /= Window
        poseX /= Window
        poseY /= Window
        poseZ /= Window
    }

    // apply corrections
    fieldX -= scan.fieldOffset.x
    fieldY = (fieldY - scan.fieldOffset.y) * scan.fieldScaleY
    fieldZ = (fieldZ - scan.fieldOffset.x) * scan.fieldScaleZ
    poseX -= poseOffset.x
    poseY = (poseY - poseOffset.y) * poseScaleY
    poseZ = (poseZ - poseOffset.z) * poseScaleY
    return new Reading(fieldX, fieldY, fieldZ, poseX, poseY, poseZ)
}

// adopt extrnally-measured calibration (for some microbits I have known...)
function characteriseAccelerometer() {
    let myName = control.deviceName()
    let dx = 0
    let dy = 0
    let dz = 0
    switch (myName) {
        case "sim-":
            poseScaleY = 1
            poseScaleZ = 1
            dx = 0
            dy = 0
            dz = 0
            break

        case "zapop":
            poseScaleY = 1042.89 / 1007.23
            poseScaleZ = 1042.89 / 992.73
            dx = -70.92
            dy = 44.597
            dz = 6.804
            break

        case "gateg":
            poseScaleY = 1017.578 / 996.736
            poseScaleZ = 1017.578 / 1026.315
            dx = -25.411
            dy = -3.251
            dz = -1.300
            break

        case "gigav":
            poseScaleY = 1057.89 / 1023.98
            poseScaleZ = 1057.89 / 1074.06
            dx = -85.33
            dy = 7.22
            dz = -18.94
            break

        case "zavov":
            poseScaleY = 1049.285 / 1059.746
            poseScaleZ = 1049.285 / 986.272
            dx = -74.082
            dy = 8.455
            dz = -7.617
            break

        default: // presume perfection until proved otherwise!
            poseScaleY = 1
            poseScaleZ = 1
            poseOffset.x = 0
            poseOffset.y = 0
            poseOffset.z = 0
            break
    } 
    poseOffset = new Vector(dx, dy, dz)
}

function simulateScan(dataset: string) {
    let times: number[]
    let samples: Vector[] = []
    let testReading: Reading
    let scanX: number[] = []
    let scanY: number[] = []
    let scanZ: number[] = []
    let testFieldX: number[] = []
    let testFieldY: number[] = []
    let testFieldZ: number[] = []
    let testPoseX: number[] = []
    let testPoseY: number[] = []
    let testPoseZ: number[] = []
    switch (dataset) {

        case "T07141743_blup70": // bottom-left upwards; dip=70
            times = [32009, 32057, 32073, 32089, 32105, 32121, 32137, 32193, 32209, 32225, 32241, 32257, 32273, 32289, 32305, 32361, 32377, 32393, 32409, 32425, 32441, 32457, 32473, 32529, 32545, 32561, 32577, 32593, 32609, 32625, 32713, 32729, 32745, 32761, 32777, 32793, 32809, 32825, 32885, 32901, 32917, 32933, 32949, 32965, 32981, 33037, 33053, 33069, 33085, 33101, 33117, 33133, 33149, 33205, 33221, 33237, 33253, 33269, 33285, 33301, 33385, 33401, 33417, 33433, 33449, 33465, 33481, 33497, 33553, 33569, 33585, 33601, 33617, 33633, 33649, 33665, 33721, 33737, 33753, 33769, 33785, 33801, 33817, 33873, 33889, 33905, 33921, 33937, 33953, 33969, 33985, 34069, 34085, 34101, 34117, 34133, 34149, 34165, 34193, 34225, 34241, 34257, 34273, 34289, 34305, 34321, 34381, 34397, 34413, 34429, 34445, 34461, 34477, 34493, 34549, 34565, 34581, 34597, 34613, 34629, 34645, 34729, 34745, 34761, 34777, 34793, 34809, 34825, 34841, 34897, 34913, 34929, 34945, 34961, 34977, 34993, 35049, 35065, 35081, 35097, 35113, 35129, 35145, 35161, 35217, 35233, 35249, 35265, 35281, 35297, 35313, 35329, 35413, 35429, 35445, 35461, 35477, 35493, 35509, 35565, 35581, 35597, 35613, 35629, 35645, 35661, 35677, 35733, 35749, 35765, 35781, 35797, 35813, 35829, 35889, 35905, 35921, 35941, 35957, 35977, 35993, 36009, 36093, 36109, 36125, 36141, 36157, 36173, 36189, 36205, 36265, 36285, 36301, 36317, 36333, 36349, 36365, 36425, 36441, 36457, 36473, 36489, 36509, 36525, 36541, 36601, 36617, 36633, 36649, 36665, 36681, 36697, 36717, 36801, 36817, 36833, 36849, 36865, 36881, 36897, 36957, 36973, 36993, 37009, 37025, 37041, 37057, 37073, 37133, 37149, 37165, 37185, 37201, 37217, 37233, 37293, 37309, 37325, 37341, 37357, 37377, 37393, 37409, 37513, 37529, 37545, 37561, 37577, 37597, 37613, 37629, 37689, 37705, 37721, 37737, 37753, 37773, 37789, 37849]
            scanX = [887.59, 889.13, 889.71, 890.29, 890.92, 891.59, 892.27, 894.78, 895.5, 896.22, 896.95, 897.69, 898.51, 899.29, 899.92, 902.1, 902.74, 903.37, 903.98, 904.58, 905.16, 905.7, 906.19, 907.75, 908.14, 908.47, 908.75, 909.03, 909.29, 909.5, 910.21, 910.23, 910.19, 910.13, 910.01, 909.81, 909.59, 909.38, 908.27, 907.95, 907.63, 907.23, 906.74, 906.24, 905.79, 904.17, 903.63, 903.04, 902.39, 901.71, 901.1, 900.47, 899.78, 897.25, 896.48, 895.71, 894.94, 894.18, 893.42, 892.63, 888.97, 888.33, 887.69, 887.06, 886.45, 885.88, 885.39, 884.91, 883.48, 883.16, 882.84, 882.57, 882.36, 882.23, 882.15, 882.08, 881.96, 881.99, 882.13, 882.31, 882.48, 882.68, 882.87, 883.82, 884.23, 884.68, 885.14, 885.63, 886.13, 886.65, 887.19, 890.34, 891.01, 891.75, 892.55, 893.35, 894.14, 894.94, 896.39, 898.01, 898.79, 899.56, 900.32, 901.05, 901.8, 902.54, 904.97, 905.53, 906.04, 906.5, 906.89, 907.31, 907.76, 908.13, 909.11, 909.33, 909.49, 909.6, 909.69, 909.74, 909.72, 909.25, 909.13, 908.98, 908.75, 908.48, 908.18, 907.89, 907.58, 906.27, 905.85, 905.39, 904.88, 904.39, 903.86, 903.32, 901.4, 900.8, 900.13, 899.43, 898.78, 898.14, 897.5, 896.82, 894.41, 893.76, 893.12, 892.48, 891.85, 891.17, 890.49, 889.84, 886.85, 886.33, 885.85, 885.4, 884.95, 884.48, 884.01, 882.74, 882.48, 882.23, 882.03, 881.87, 881.7, 881.53, 881.44, 881.36, 881.36, 881.42, 881.55, 881.72, 881.89, 882.06, 883.09, 883.45, 883.83, 884.43, 884.98, 885.72, 886.35, 886.94, 890.46, 891.23, 891.97, 892.68, 893.45, 894.29, 895.16, 896.04, 899.26, 900.29, 901.1, 901.88, 902.62, 903.3, 903.95, 906.16, 906.66, 907.11, 907.53, 907.91, 908.32, 908.62, 908.92, 909.63, 909.73, 909.74, 909.7, 909.68, 909.67, 909.62, 909.46, 908.19, 907.91, 907.64, 907.3, 906.94, 906.64, 906.29, 904.66, 904.21, 903.61, 903.11, 902.59, 902.04, 901.47, 900.89, 898.83, 898.25, 897.66, 896.93, 896.36, 895.78, 895.17, 892.83, 892.24, 891.69, 891.13, 890.52, 889.72, 889.11, 888.53, 885.11, 884.67, 884.28, 883.88, 883.49, 883.08, 882.74, 882.42, 881.65, 881.51, 881.38, 881.31, 881.27, 881.28, 881.32, 881.75, 881.16]
            scanY = [1586.86, 1587.98, 1588.36, 1588.68, 1589, 1589.33, 1589.65, 1590.42, 1590.58, 1590.69, 1590.77, 1590.85, 1590.91, 1590.88, 1590.8, 1590.45, 1590.3, 1590.09, 1589.86, 1589.56, 1589.22, 1588.9, 1588.58, 1587.2, 1586.76, 1586.24, 1585.67, 1585.17, 1584.69, 1584.12, 1580.57, 1579.91, 1579.24, 1578.54, 1577.81, 1577.08, 1576.37, 1575.7, 1573.32, 1572.64, 1571.95, 1571.27, 1570.54, 1569.85, 1569.19, 1567.16, 1566.65, 1566.15, 1565.65, 1565.18, 1564.74, 1564.34, 1563.96, 1562.9, 1562.71, 1562.58, 1562.46, 1562.37, 1562.32, 1562.28, 1562.83, 1563.04, 1563.27, 1563.55, 1563.85, 1564.18, 1564.57, 1564.97, 1566.67, 1567.26, 1567.86, 1568.43, 1569.04, 1569.67, 1570.34, 1571.08, 1573.56, 1574.22, 1574.9, 1575.69, 1576.44, 1577.15, 1577.85, 1580.1, 1580.72, 1581.35, 1581.96, 1582.56, 1583.19, 1583.75, 1584.3, 1586.89, 1587.3, 1587.66, 1587.97, 1588.28, 1588.58, 1588.85, 1589.2, 1589.44, 1589.49, 1589.44, 1589.33, 1589.18, 1588.97, 1588.75, 1587.81, 1587.52, 1587.15, 1586.71, 1586.19, 1585.66, 1585.17, 1584.7, 1582.82, 1582.22, 1581.61, 1581.01, 1580.42, 1579.8, 1579.15, 1575.98, 1575.36, 1574.75, 1574.14, 1573.54, 1572.96, 1572.38, 1571.81, 1569.68, 1569.09, 1568.52, 1567.99, 1567.52, 1567.02, 1566.5, 1564.87, 1564.47, 1564.16, 1563.89, 1563.54, 1563.15, 1562.83, 1562.58, 1562.06, 1561.96, 1561.89, 1561.84, 1561.78, 1561.75, 1561.81, 1561.93, 1562.98, 1563.29, 1563.62, 1563.94, 1564.28, 1564.64, 1565.02, 1566.62, 1567.08, 1567.55, 1568.05, 1568.58, 1569.14, 1569.68, 1570.18, 1572.16, 1572.79, 1573.48, 1574.17, 1574.85, 1575.52, 1576.18, 1578.69, 1579.38, 1580.11, 1581.01, 1581.7, 1582.52, 1583.15, 1583.72, 1586.25, 1586.71, 1587.18, 1587.59, 1587.97, 1588.29, 1588.54, 1588.75, 1589.03, 1588.98, 1588.89, 1588.76, 1588.61, 1588.43, 1588.2, 1586.94, 1586.51, 1586.08, 1585.67, 1585.27, 1584.71, 1584.17, 1583.53, 1581.21, 1580.6, 1579.94, 1579.27, 1578.61, 1577.96, 1577.31, 1576.5, 1573.1, 1572.49, 1571.93, 1571.34, 1570.65, 1569.96, 1569.34, 1567.43, 1566.95, 1566.37, 1565.97, 1565.57, 1565.15, 1564.73, 1564.34, 1563.05, 1562.75, 1562.5, 1562.25, 1562.04, 1561.85, 1561.7, 1561.32, 1561.29, 1561.28, 1561.3, 1561.35, 1561.48, 1561.58, 1561.71, 1563.58, 1564, 1564.44, 1564.89, 1565.37, 1565.9, 1566.32, 1566.79, 1568.82, 1569.37, 1569.94, 1570.58, 1571.21, 1571.97, 1572.6, 1575.1, 1566.09]
            scanZ = [424.65, 424.91, 425.05, 425.15, 425.24, 425.37, 425.57, 426.47, 426.72, 426.93, 427.14, 427.37, 427.62, 427.9, 428.19, 429.33, 429.66, 429.95, 430.25, 430.57, 430.88, 431.19, 431.54, 432.76, 433.13, 433.54, 433.92, 434.26, 434.56, 434.86, 436.42, 436.7, 436.97, 437.25, 437.52, 437.77, 437.99, 438.19, 438.85, 439.02, 439.16, 439.28, 439.38, 439.46, 439.57, 439.72, 439.7, 439.67, 439.6, 439.47, 439.38, 439.35, 439.31, 438.83, 438.62, 438.39, 438.15, 437.87, 437.58, 437.27, 435.39, 435.02, 434.7, 434.36, 434, 433.62, 433.21, 432.78, 431.34, 430.96, 430.6, 430.24, 429.84, 429.48, 429.16, 428.84, 427.79, 427.51, 427.26, 427, 426.75, 426.57, 426.44, 426.01, 425.91, 425.82, 425.7, 425.61, 425.55, 425.49, 425.45, 425.69, 425.8, 425.95, 426.12, 426.33, 426.55, 426.74, 427.16, 427.75, 428.05, 428.39, 428.81, 429.2, 429.55, 429.91, 431.31, 431.67, 432.02, 432.38, 432.72, 433.09, 433.45, 433.81, 435.05, 435.36, 435.68, 436.02, 436.32, 436.61, 436.92, 438.33, 438.53, 438.77, 438.97, 439.11, 439.29, 439.48, 439.61, 439.84, 439.93, 439.97, 439.96, 439.96, 439.95, 439.97, 439.8, 439.67, 439.55, 439.47, 439.33, 439.2, 439.08, 438.92, 438.25, 438.04, 437.84, 437.57, 437.28, 437, 436.68, 436.35, 434.66, 434.3, 433.91, 433.56, 433.23, 432.91, 432.58, 431.48, 431.21, 430.9, 430.59, 430.26, 429.91, 429.61, 429.34, 428.46, 428.2, 427.98, 427.74, 427.5, 427.33, 427.16, 426.47, 426.32, 426.18, 426, 425.91, 425.8, 425.74, 425.75, 426.1, 426.22, 426.34, 426.47, 426.66, 426.88, 427.1, 427.35, 428.54, 428.98, 429.35, 429.7, 430.04, 430.41, 430.79, 432.22, 432.61, 432.99, 433.35, 433.71, 434.16, 434.5, 434.81, 436.04, 436.39, 436.7, 436.95, 437.26, 437.6, 437.9, 438.22, 439.05, 439.15, 439.24, 439.37, 439.52, 439.68, 439.78, 439.89, 439.95, 440.02, 440.02, 439.95, 439.88, 439.85, 439.85, 439.58, 439.43, 439.29, 439.12, 438.97, 438.8, 438.59, 437.8, 437.55, 437.27, 437.06, 436.82, 436.41, 436.1, 435.79, 433.47, 433.13, 432.81, 432.47, 432.13, 431.72, 431.38, 431.05, 429.9, 429.59, 429.32, 429.04, 428.72, 428.38, 428.13, 427.27, 430.74]
            testFieldX = [881.04, 880.44, 889.41, 901.18, 910.09, 911.06, 901.67, 889.44, 880.74, 880.39, 888.66, 900.99, 910.05, 910.09, 901.37, 889.26, 880.29, 879.88, 888.69, 900.51, 909.99, 909.77, 901.22, 888.58, 879.79]
            testFieldY = [1566.06, 1577.64, 1588.18, 1591.76, 1585.86, 1573.95, 1562.89, 1559.31, 1565.21, 1576.89, 1587.86, 1591.16, 1584.79, 1573.18, 1562.72, 1559.14, 1565.25, 1576.5, 1587.41, 1590.28, 1584.86, 1572.92, 1562.46, 1558.41, 1564.11]
            testFieldZ = [430.54, 425.51, 424.33, 428.21, 434.21, 439.48, 440.53, 437.04, 430.59, 425.85, 424.91, 428.08, 434.46, 439.11, 440.31, 436.91, 430.22, 425.46, 424.44, 427.63, 434.36, 439.48, 440.25, 436.84, 430.5]
            // poses were never captured...
            testPoseX = []
            testPoseY = []
            testPoseZ = []    
            break

        case "T07260757_dash70": // angled forward like a dash-board: dip=70
            times = [9229, 9245, 9261, 9277, 9293, 9309, 9325, 9341, 9357, 9373, 9389, 9405, 9421, 9437, 9453, 9469, 9485, 9501, 9517, 9533, 9549, 9565, 9581, 9597, 9613, 9629, 9645, 9661, 9677, 9693, 9709, 9725, 9741, 9757, 9773, 9789, 9805, 9821, 9837, 9853, 9869, 9885, 9901, 9917, 9933, 9949, 9965, 9981, 9997, 10013, 10029, 10045, 10061, 10077, 10093, 10109, 10125, 10141, 10157, 10173, 10189, 10205, 10221, 10237, 10253, 10269, 10285, 10301, 10317, 10333, 10349, 10365, 10381, 10397, 10413, 10429, 10445, 10461, 10477, 10493, 10509, 10525, 10541, 10557, 10573, 10589, 10605, 10621, 10637, 10653, 10669, 10685, 10701, 10717, 10733, 10749, 10765, 10781, 10797, 10813, 10829, 10845, 10861, 10877, 10893, 10909, 10925, 10941, 10957, 10973, 10989, 11005, 11021, 11037, 11053, 11069, 11085, 11101, 11117, 11133, 11149, 11165, 11181, 11197, 11213, 11229, 11245, 11261, 11277, 11293, 11309, 11325, 11341, 11357, 11373, 11389, 11405, 11421, 11437, 11453, 11469, 11485, 11501, 11517, 11533, 11549, 11565, 11581, 11597, 11613, 11629, 11645, 11661, 11677, 11693, 11709, 11725, 11741, 11757, 11773, 11789, 11805, 11821, 11837, 11853, 11869, 11885, 11901, 11917, 11933, 11949, 11965, 11981, 11997, 12013, 12029, 12045, 12061, 12077, 12093, 12109, 12125, 12141, 12157, 12173, 12189, 12205, 12221, 12237, 12253, 12269, 12285, 12301, 12317, 12333, 12349, 12365, 12381, 12397, 12413, 12429, 12445, 12461, 12477, 12493, 12509, 12525, 12541, 12557, 12573, 12589, 12605, 12621, 12637, 12653, 12669, 12685, 12701, 12717, 12733, 12749, 12765, 12781, 12797, 12813, 12829, 12845, 12861, 12877, 12893, 12909, 12925, 12941, 12957, 12973, 12989, 13005, 13021, 13037, 13053, 13069, 13085, 13101, 13117, 13133, 13149, 13165, 13181, 13197, 13213, 13229, 13245, 13261, 13277, 13293, 13309, 13325, 13341, 13357, 13373, 13389, 13405, 13421, 13437, 13453, 13469, 13485, 13501, 13517, 13533, 13549, 13565, 13581, 13597, 13613, 13629, 13645, 13661, 13677, 13693, 13709, 13725, 13741, 13757, 13773, 13789, 13805, 13821, 13837, 13853, 13869, 13885, 13901, 13917, 13933, 13949, 13965, 13981, 13997, 14013, 14029, 14045, 14061, 14077, 14093, 14109, 14125, 14141, 14157, 14173, 14189, 14205, 14221, 14237, 14253, 14269, 14285, 14301, 14317, 14333, 14349, 14365, 14381, 14397, 14413, 14429, 14445, 14461, 14477, 14493, 14509, 14525, 14541, 14557, 14573, 14589, 14605, 14621, 14637, 14653, 14669, 14685, 14701, 14717, 14733, 14749, 14765, 14781, 14797, 14813, 14829, 14845, 14861, 14877, 14893, 14909, 14925, 14941, 14957, 14973, 14989, 15005, 15021, 15037]
            scanX = [-17.069, -17.374, -17.698, -17.999, -18.321, -18.653, -18.97, -19.305, -19.637, -19.946, -20.281, -20.632, -20.946, -21.233, -21.501, -21.799, -22.115, -22.405, -22.722, -23.01, -23.227, -23.506, -23.784, -24.013, -24.254, -24.481, -24.675, -24.865, -25.065, -25.24, -25.394, -25.534, -25.669, -25.805, -25.949, -26.127, -26.288, -26.402, -26.507, -26.576, -26.6, -26.635, -26.693, -26.753, -26.807, -26.897, -27.024, -27.14, -27.245, -27.315, -27.347, -27.375, -27.408, -27.445, -27.484, -27.507, -27.509, -27.544, -27.581, -27.543, -27.469, -27.455, -27.441, -27.383, -27.349, -27.332, -27.309, -27.276, -27.199, -27.07, -26.919, -26.769, -26.638, -26.526, -26.377, -26.215, -26.094, -25.989, -25.848, -25.64, -25.398, -25.171, -24.977, -24.787, -24.56, -24.268, -24.003, -23.779, -23.55, -23.32, -23.084, -22.81, -22.521, -22.203, -21.849, -21.559, -21.254, -20.951, -20.632, -20.258, -19.95, -19.635, -19.215, -18.802, -18.43, -18.065, -17.698, -17.325, -16.951, -16.557, -16.159, -15.782, -15.359, -14.895, -14.472, -14.052, -13.637, -13.211, -12.771, -12.323, -11.882, -11.453, -11.029, -10.609, -10.146, -9.633, -9.085, -8.594, -8.184, -7.743, -7.243, -6.761, -6.341, -5.97, -5.545, -5.107, -4.74, -4.364, -3.99, -3.654, -3.301, -2.941, -2.612, -2.312, -2.057, -1.838, -1.618, -1.392, -1.145, -0.874, -0.615, -0.397, -0.212, -0.036, 0.158, 0.32, 0.423, 0.511, 0.631, 0.772, 0.881, 0.975, 1.06, 1.144, 1.208, 1.233, 1.283, 1.368, 1.405, 1.44, 1.458, 1.382, 1.281, 1.195, 1.11, 1.053, 0.98, 0.862, 0.747, 0.642, 0.494, 0.3, 0.116, -0.084, -0.3, -0.536, -0.778, -0.985, -1.174, -1.4, -1.658, -1.92, -2.231, -2.58, -2.905, -3.182, -3.428, -3.706, -4.038, -4.38, -4.752, -5.136, -5.52, -5.942, -6.387, -6.792, -7.182, -7.577, -7.957, -8.337, -8.743, -9.178, -9.588, -9.956, -10.326, -10.748, -11.205, -11.655, -12.078, -12.491, -12.9, -13.354, -13.851, -14.339, -14.803, -15.24, -15.658, -16.086, -16.529, -17.005, -17.499, -17.915, -18.29, -18.651, -18.988, -19.308, -19.653, -20, -20.33, -20.71, -21.109, -21.48, -21.791, -22.063, -22.324, -22.609, -22.921, -23.212, -23.504, -23.783, -24.031, -24.291, -24.554, -24.767, -24.963, -25.193, -25.457, -25.656, -25.773, -25.894, -26.043, -26.219, -26.385, -26.517, -26.625, -26.71, -26.788, -26.859, -26.917, -26.977, -27.018, -27.043, -27.038, -26.976, -26.879, -26.785, -26.679, -26.57, -26.489, -26.396, -26.238, -26.056, -25.883, -25.678, -25.42, -25.152, -24.907, -24.61, -24.295, -23.987, -23.634, -23.272, -22.93, -22.61, -22.272, -21.845, -21.398, -20.943, -20.463, -20.005, -19.561, -19.116, -18.585, -18.01, -17.518, -17.101, -16.654, -16.169, -15.653, -15.124, -14.608, -14.118, -13.644, -13.1, -12.565, -12.11, -11.654, -11.146, -10.617, -10.109, -9.632, -9.163, -8.65, -8.137, -7.649, -7.164, -6.708, -6.279, -5.835, -5.397, -4.977, -4.522, -4.062, -3.654, -3.276, -2.938, -2.588, -2.205, -1.882, -1.606, -1.321, -1.028, -0.748, -0.486, -0.254, -0.048, 0.159, 0.382, 0.597, 0.748, 0.859, 0.968, 1.073, 1.163, 1.218, 1.243, 1.179, 1.149, 1.193, 1.206, 1.206, 1.179, 1.122, 1.036]
            scanY = [-4.611, -4.555, -4.504, -4.438, -4.365, -4.312, -4.302, -4.26, -4.176, -4.106, -3.995, -3.876, -3.76, -3.595, -3.425, -3.283, -3.16, -3.048, -2.932, -2.78, -2.616, -2.466, -2.322, -2.152, -1.957, -1.773, -1.625, -1.516, -1.37, -1.218, -1.131, -1.072, -0.957, -0.815, -0.698, -0.551, -0.425, -0.349, -0.232, -0.081, 0.044, 0.198, 0.424, 0.662, 0.826, 0.944, 1.065, 1.215, 1.397, 1.543, 1.704, 1.929, 2.147, 2.326, 2.49, 2.653, 2.847, 3.051, 3.236, 3.452, 3.688, 3.95, 4.2, 4.404, 4.635, 4.886, 5.127, 5.348, 5.556, 5.785, 5.993, 6.164, 6.36, 6.598, 6.801, 7.007, 7.26, 7.503, 7.719, 7.93, 8.155, 8.399, 8.648, 8.876, 9.056, 9.256, 9.515, 9.738, 9.91, 10.093, 10.306, 10.466, 10.607, 10.801, 10.998, 11.177, 11.367, 11.554, 11.737, 11.93, 12.067, 12.193, 12.36, 12.505, 12.651, 12.817, 12.982, 13.138, 13.242, 13.33, 13.437, 13.534, 13.596, 13.634, 13.641, 13.649, 13.688, 13.724, 13.784, 13.844, 13.873, 13.919, 13.921, 13.858, 13.787, 13.763, 13.772, 13.742, 13.694, 13.605, 13.497, 13.406, 13.283, 13.144, 13.025, 12.895, 12.792, 12.649, 12.442, 12.259, 12.059, 11.837, 11.642, 11.417, 11.14, 10.874, 10.645, 10.416, 10.148, 9.884, 9.649, 9.455, 9.244, 8.973, 8.714, 8.443, 8.151, 7.904, 7.659, 7.363, 7.067, 6.817, 6.601, 6.351, 6.058, 5.791, 5.527, 5.254, 4.991, 4.712, 4.458, 4.199, 3.924, 3.686, 3.469, 3.237, 3, 2.784, 2.533, 2.26, 2.01, 1.762, 1.489, 1.195, 0.907, 0.651, 0.424, 0.21, -0.005, -0.175, -0.29, -0.443, -0.653, -0.911, -1.159, -1.373, -1.58, -1.73, -1.861, -2.053, -2.22, -2.328, -2.492, -2.707, -2.91, -3.088, -3.217, -3.316, -3.44, -3.578, -3.685, -3.776, -3.925, -4.12, -4.257, -4.33, -4.378, -4.402, -4.41, -4.438, -4.483, -4.524, -4.587, -4.645, -4.672, -4.663, -4.662, -4.667, -4.653, -4.646, -4.633, -4.601, -4.547, -4.483, -4.399, -4.303, -4.242, -4.206, -4.126, -4.001, -3.868, -3.74, -3.603, -3.451, -3.301, -3.15, -2.966, -2.759, -2.597, -2.428, -2.223, -2.008, -1.794, -1.571, -1.329, -1.055, -0.819, -0.616, -0.392, -0.154, 0.122, 0.413, 0.711, 1.036, 1.365, 1.692, 1.996, 2.314, 2.662, 2.98, 3.297, 3.638, 3.989, 4.289, 4.563, 4.846, 5.157, 5.493, 5.806, 6.105, 6.424, 6.754, 7.055, 7.35, 7.671, 8.01, 8.333, 8.647, 8.913, 9.184, 9.482, 9.735, 9.959, 10.177, 10.442, 10.702, 10.94, 11.164, 11.342, 11.55, 11.793, 12.047, 12.279, 12.464, 12.673, 12.896, 13.026, 13.101, 13.201, 13.306, 13.408, 13.521, 13.584, 13.61, 13.693, 13.807, 13.885, 13.903, 13.887, 13.867, 13.848, 13.777, 13.661, 13.555, 13.452, 13.312, 13.145, 13.012, 12.877, 12.715, 12.528, 12.304, 12.081, 11.898, 11.703, 11.453, 11.239, 11.083, 10.887, 10.657, 10.44, 10.209, 9.912, 9.583, 9.27, 8.965, 8.639, 8.317, 7.99, 7.645, 7.345, 7.076, 6.785, 6.494, 6.241, 5.944, 5.608, 5.281, 5.004, 4.743, 4.46, 4.175, 3.862, 3.524]
            scanZ = [80.422, 80.44, 80.463, 80.47, 80.444, 80.361, 80.313, 80.314, 80.308, 80.25, 80.165, 80.1, 80.028, 79.982, 79.966, 79.958, 79.938, 79.913, 79.89, 79.856, 79.825, 79.763, 79.728, 79.707, 79.663, 79.636, 79.571, 79.493, 79.414, 79.357, 79.317, 79.254, 79.197, 79.191, 79.172, 79.115, 79.11, 79.131, 79.13, 79.123, 79.11, 79.045, 78.954, 78.871, 78.795, 78.758, 78.704, 78.63, 78.588, 78.57, 78.519, 78.408, 78.321, 78.233, 78.121, 78.052, 77.988, 77.897, 77.788, 77.676, 77.586, 77.478, 77.346, 77.2, 77.087, 77.022, 76.997, 76.939, 76.882, 76.835, 76.74, 76.654, 76.561, 76.478, 76.374, 76.301, 76.239, 76.131, 76.045, 75.96, 75.863, 75.768, 75.675, 75.583, 75.519, 75.447, 75.329, 75.242, 75.195, 75.105, 75.051, 75.033, 74.957, 74.883, 74.827, 74.764, 74.693, 74.588, 74.533, 74.561, 74.506, 74.392, 74.308, 74.178, 74.04, 73.954, 73.921, 73.89, 73.764, 73.643, 73.585, 73.538, 73.469, 73.433, 73.432, 73.385, 73.317, 73.32, 73.335, 73.297, 73.283, 73.322, 73.356, 73.356, 73.342, 73.337, 73.331, 73.306, 73.284, 73.317, 73.398, 73.446, 73.444, 73.423, 73.416, 73.47, 73.571, 73.663, 73.727, 73.773, 73.861, 73.967, 74.097, 74.241, 74.374, 74.52, 74.66, 74.742, 74.804, 74.873, 74.911, 74.965, 75.086, 75.283, 75.454, 75.585, 75.717, 75.851, 75.94, 76.016, 76.144, 76.253, 76.327, 76.414, 76.48, 76.548, 76.709, 76.913, 77.053, 77.143, 77.217, 77.279, 77.389, 77.534, 77.633, 77.711, 77.842, 77.98, 78.053, 78.088, 78.126, 78.194, 78.302, 78.42, 78.527, 78.592, 78.599, 78.639, 78.681, 78.734, 78.84, 78.924, 78.974, 79.013, 79.076, 79.072, 78.992, 79.025, 79.105, 79.126, 79.186, 79.279, 79.337, 79.384, 79.467, 79.555, 79.629, 79.696, 79.751, 79.808, 79.826, 79.857, 79.944, 80.023, 80.082, 80.07, 80.034, 80.033, 80.039, 80.044, 80.076, 80.177, 80.214, 80.192, 80.21, 80.253, 80.305, 80.353, 80.43, 80.465, 80.447, 80.463, 80.456, 80.399, 80.38, 80.395, 80.428, 80.392, 80.308, 80.283, 80.267, 80.239, 80.208, 80.205, 80.195, 80.141, 80.098, 80.072, 80.055, 79.993, 79.951, 79.897, 79.771, 79.669, 79.608, 79.533, 79.405, 79.302, 79.241, 79.205, 79.132, 79.026, 78.93, 78.816, 78.705, 78.571, 78.443, 78.325, 78.179, 78.032, 77.891, 77.723, 77.571, 77.498, 77.416, 77.297, 77.175, 76.997, 76.823, 76.71, 76.593, 76.467, 76.398, 76.321, 76.179, 76.064, 75.984, 75.871, 75.749, 75.691, 75.61, 75.471, 75.338, 75.191, 75.036, 74.889, 74.801, 74.773, 74.74, 74.664, 74.536, 74.406, 74.281, 74.16, 74.057, 73.989, 73.919, 73.85, 73.813, 73.786, 73.741, 73.681, 73.679, 73.655, 73.596, 73.575, 73.549, 73.498, 73.473, 73.482, 73.483, 73.505, 73.524, 73.555, 73.576, 73.609, 73.691, 73.806, 73.943, 74.042, 74.099, 74.105, 74.079, 74.107, 74.212, 74.347, 74.449, 74.51, 74.619, 74.771, 74.869, 74.953, 75.108, 75.295, 75.446, 75.531, 75.577, 75.671, 75.796, 75.897, 76.082, 76.261, 76.348, 76.42, 76.479, 76.587, 76.723, 76.868, 77.032, 77.163, 77.262, 77.379, 77.512, 77.613]
            testFieldX = [1.093, 1.029, 1.05, 1.929, 2.014, 2.164, 2.014, 2.229, 2.079, 0.814, 1.071, 1.05, -1.629, -1.157, -1.543, -4.221, -4.221, -4.2, -7.307, -7.136, -7.136, -10.543, -10.821, -10.843, -14.679, -14.743, -14.464, -17.85, -18.279, -17.743, -21.236, -21.15, -21.129, -24, -23.871, -24.107, -26.271, -26.164, -26.164, -27.45, -27.471, -27.536, -28.2, -27.879, -27.857, -27.3, -27.129, -27.686, -25.714, -25.479, -25.629, -22.714, -22.564, -22.8, -19.243, -19.5, -19.457, -15.321, -15.557, -14.914, -10.8, -10.929, -11.25, -7.5, -7.564, -7.436, -3.664, -3.75, -3.857, -1.071, -1.071, -0.621, 0.943, 0.814, 0.643]
            testFieldY = [9.086, 8.764, 8.764, 6.664, 6.664, 6.643, 3.621, 3.964, 4.05, 1.671, 1.757, 1.543, -0.321, -0.793, -0.579, -2.229, -2.636, -2.636, -3.621, -3.6, -3.836, -4.5, -4.414, -4.457, -4.993, -5.143, -4.8, -4.607, -5.057, -4.714, -3.9, -3.9, -3.471, -2.186, -2.55, -2.486, -0.579, -0.514, -0.621, 1.779, 1.779, 1.671, 4.05, 4.114, 3.857, 6.707, 6.193, 6.321, 8.721, 9.3, 8.871, 11.421, 11.336, 11.336, 12.836, 12.921, 12.686, 14.4, 14.486, 14.014, 14.679, 14.914, 14.55, 14.379, 14.293, 14.336, 12.729, 12.879, 13.2, 11.271, 10.864, 10.671, 9.043, 8.743, 8.7]
            testFieldZ = [75.514, 75.557, 75.043, 75.943, 76.093, 76.071, 77.293, 77.314, 77.271, 78.343, 78.236, 77.914, 78.771, 78.557, 78.771, 78.986, 79.114, 79.286, 79.35, 79.693, 79.371, 79.886, 79.779, 79.95, 80.207, 79.8, 79.971, 80.079, 79.864, 80.164, 79.929, 80.186, 80.464, 79.907, 79.929, 79.843, 79.65, 79.221, 79.393, 78.236, 78.214, 78.45, 77.336, 77.379, 77.55, 76.157, 76.35, 76.564, 75.364, 75.043, 75.664, 74.4, 74.229, 74.379, 73.95, 73.95, 73.629, 73.264, 73.05, 73.221, 72.879, 72.6, 73.157, 72.514, 72.643, 72.964, 73.779, 73.5, 73.457, 74.186, 74.529, 74.529, 74.979, 75.086, 75]
            // poses were never captured...
            testPoseX = []
            testPoseY = []
            testPoseZ = []
            break
    }
    
    // transpose the three arrays into the scan array of triples...
    scan.samples = []
    scan.times = []
    for (let i = 0; i < times.length; i++) {
        scan.samples.push(new Vector(scanX[i], scanY[i], scanZ[i]))
        scan.times.push(times[i])
    }

    // assemble the array of test readings...
    testReadings = []
    for (let n = 0; n < testFieldX.length; n++) {
        if (testPoseX.length > 0 ) {
            testReading = new Reading(
                testFieldX[n], testFieldY[n], testFieldZ[n], 
                testPoseX[n], testPoseY[n], testPoseZ[n])
        } else { // this is an old test dataset, for which pose data was never captured
            testReading = new Reading(testFieldX[n], testFieldY[n], testFieldZ[n], 0, 0, -1023) // for now, always pretend it was face-up!
        }
        testReadings.push(testReading)
    }
    testIndex = 0
    return 0 // never fails!
}



// =============== FOREGROUND CODE =================
scan = new Scan()
basic.clearScreen()
basic.showString(control.deviceName())

/*
* We are using three different 3D frames of reference:
*
*       XYZ: the microbit Sensor-Frame
*       RFD: the buggy Body-Frame (Right, Front, Down)
*       ENG: the World-Frame in which it is navigating (East, North, Gravity)
*
*/
// await button-pressing...

