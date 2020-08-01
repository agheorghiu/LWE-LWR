namespace RSPRounding {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics; 
    open Microsoft.Quantum.Convert;

    // controlled version of the previous function (controlled on b)
    operation controlledAddConstant(b : Qubit, y : Int, modulus : Int, res : Qubit[]) : Unit {
        body (...) {
            let adder = addConstantToRegister(y, modulus, _);
            (Controlled adder)([b], res);
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;        
    }

    // coherently adding <c, x> + b . y mod modulus to the result register
    operation coherentLWE(b : Qubit, x : (Qubit[])[], res : Qubit[], y : Int, c : Int[], modulus : Int) : Unit {
        body (...) {
            let n = Length(c);
            for (j in 0 .. n - 1) {
                if (c[j] != 0) {
                    addRegisters(x[j], res, c[j], modulus); // |x>|cx + res>
                }
            }
            controlledAddConstant(b, y, modulus, res); // |b>|x>|cx + by + res>
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;
    }

    // take most significant bit as rounded value
    operation roundMSB(anc: Qubit[], res: Qubit) : Unit {
        let bits = Length(anc);
        CNOT(anc[0], res);
    }

    operation RSPMain (A : Int[][], y : Int[], modulus : Int, rnd : Int) : (Int, Int[], Int[]) {
        let m = Length(A);
        let n = Length(A[0]);
        let bits = Ceiling(Lg(IntAsDouble(modulus)));
        let nQubits = bits * n + bits + m;
        mutable roundRes = new Int[m];
        mutable bMeas = 0;
        mutable d = new Int[bits * n];

        using (qs = Qubit[nQubits]) {
            let b = qs[0];

            // initialize x as vector of blocks of qubits
            mutable x = new (Qubit[])[n];
            for (i in 0 .. n - 1) {
                set x w/= i <- (qs[(1 + bits * i) .. (bits * (i + 1))]);
            }

            // prep superpositions over b and x
            // have \sum_{b, x} |b>|x>|0>
            H(b);
            for (i in 0 .. n - 1) {
                ApplyToEach(H, x[i]);
            }

            // ancilla register for partial results
            let anc = qs[1 + n * bits .. (n + 1) * bits];
            ApplyToEach(H, anc);

            // register for rounded qubits
            let r = qs[(n + 1) * bits + 1 .. bits * n + bits + m - 1] + [anc[0]]; // last rounded qubit is MSB of ancilla

            // compute rounded result; use ancilla to store pre-rounded result; clear ancilla each time
            // final state should be |b>|x>|0>|round(Ax+b . (As + e'))>
            for (i in 0 .. m - 2) {
                coherentLWE(b, x, anc, y[i], A[i], modulus);
                (Adjoint QFTNormal)(anc);
                roundMSB(anc, r[i]);
                QFTNormal(anc);
                (Adjoint coherentLWE)(b, x, anc, y[i], A[i], modulus);
            }
            // in the last iteration MSB of ancilla is the rounded result
            coherentLWE(b, x, anc, y[m - 1], A[m - 1], modulus);
            (Adjoint QFTNormal)(anc);            


            // Hadamard pre-image register
            H(b);
            for (i in 0 .. n - 1) {
                ApplyToEach(H, x[i]);
            }

            // measure everything
            set bMeas = MeasureInteger(LittleEndian([b]));

            for (i in 0 .. n * bits - 1) {
                set d w/= i <- MeasureInteger(LittleEndian([qs[i + 1]]));
            }

            for (i in 0 .. m - 1) {
                set roundRes w/= i <- MeasureInteger(LittleEndian([r[i]]));
            }

            // reset qubits
            ResetAll(qs);
        }

        return (bMeas, d, roundRes);
    }
}