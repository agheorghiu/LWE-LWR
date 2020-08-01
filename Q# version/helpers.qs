namespace RSPRounding {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics; 
    open Microsoft.Quantum.Convert;

    // normal depth QFT, assumes qs is LittleEndian
    // credit to Qiskit documentation for their implementation, of which this is the iterative version
    operation QFTNormal(qs : Qubit[]) : Unit is Adj + Ctl {
        body (...) {
            H(qs[0]);
            Controlled S([qs[1]], qs[0]);
            H(qs[1]);
        }
        // body (...) {
        //     let n = Length(qs);
        //     for (i in 0 .. n - 1) {
        //         H((qs)[n - i - 1]);
        //         if (i != n - 1) {
        //             for (j in 0 .. n - i - 2) {
        //                 Controlled R1Frac([(qs)[j]], (1, n - 1 - i - j, (qs)[n - i - 1]));
        //             }
        //         }
        //     }
        //     for (i in 0 .. (n / 2) - 1) {
        //         // SWAP((qs!)[i], (qs!)[n - 1 - i]);
        //         CNOT(qs[i], qs[n - 1 - i]);
        //         CNOT(qs[n - 1 - i], (qs)[i]);
        //         CNOT(qs[i], qs[n - 1 - i]);
        //     }
        // }

        // controlled (controls, ...) {
        //     let n = Length(qs);
        //     for (i in 0 .. n - 1) {
        //         Controlled H(controls, qs[n - i - 1]);
        //         for (j in 0 .. n - i - 2) {
        //             // doubly controlled rotation
        //             Controlled R1Frac(controls + [qs[j]], (1, n - 1 - i - j, qs[n - i - 1]));
        //         }
        //     }
        //     for (i in 0 .. n / 2 - 1) {
        //         // Controlled SWAP(controls, ((qs!)[i], (qs!)[n - 1 - i]));
        //         // three Toffolis, see Beauregard paper to improve
        //         Controlled X(controls + [qs[i]], qs[n - 1 - i]);
        //         Controlled X(controls + [qs[n - 1 - i]], qs[i]);
        //         Controlled X(controls + [qs[i]], qs[n - 1 - i]);
        //     }
        // }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    // maps |x>|y> to |cx + y (mod modulus)>
    // this version only works if modulus = 2^n, where n = Length(x) = Length(y)
    operation addRegisters(x : Qubit[], y : Qubit[], c : Int, modulus : Int) : Unit is Adj + Ctl{
        let n = Length(x);
        //QFTNormal(y);
        for (i in 0 .. n - 1) {
            let sum = ((2 ^ i) * c) % (modulus);
            for (j in 0 .. n - 1) {
                if (sum % (2 ^ (n - 1 - j)) != 0 or (sum % (2 ^ (n - 1 - j)) == 0 and (sum / (2 ^ (n - 1 - j))) % 2 != 0)) {
                    Controlled R1Frac([x[i]], (sum, (n - 1) - j, y[j]));
                }
            }
        }
        //Adjoint QFTNormal(y);
    }

    // maps |x> to |x + c (mod modulus)>
    // this version only works if modulus = 2^n, where n = Length(x) = Length(y)
    operation addConstantToRegister(c : Int, modulus : Int, res : Qubit[]) : Unit is Adj + Ctl {
        body (...) {
            let n = Length(res);
            //QFTNormal(res);
            for (i in 0 .. n - 1) {
                if (c % (2 ^ (n - 1 - i)) != 0 or (c % (2 ^ (n - 1 - i)) == 0 and (c / (2 ^ (n - 1 - i))) % 2 != 0)) {
                    R1Frac(c, (n - 1) - i, (res)[i]);
                }
            }
            //Adjoint QFTNormal(res);
        }

        controlled (controls, ...) {
            let n = Length(res);
            //Controlled QFTNormal(controls, res);
            for (i in 0 .. n - 1) {
                if (c % (2 ^ (n - 1 - i)) != 0 or (c % (2 ^ (n - 1 - i)) == 0 and (c / (2 ^ (n - 1 - i))) % 2 != 0)) {
                    Controlled R1Frac(controls, (c, (n - 1) - i, (res)[i]));
                }
            }
            //Controlled (Adjoint QFTNormal)(controls, res);
        }
        adjoint invert;
    }
}