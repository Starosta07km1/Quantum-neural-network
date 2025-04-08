import Std.Diagnostics.*;

operation Main() : (Result, Result)[] {
    let bellStateTuples = [
        ("|Φ+〉", PreparePhiPlus),
        ("|Φ-〉", PreparePhiMinus),
        ("|Ψ+〉", PreparePsiPlus),
        ("|Ψ-〉", PreparePsiMinus)
    ];

    mutable measurements = [];
    for (label, prepare) in bellStateTuples {
        use register = Qubit[2];
        prepare(register);
        Message($"Bell state {label}:");
        DumpMachine();
        set measurements += [(MResetZ(register[0]), MResetZ(register[1]))];
    }
    return measurements;
}

operation PreparePhiPlus(register : Qubit[]) : Unit {
    H(register[0]);                 // |+0〉
    CNOT(register[0], register[1]); // 1/sqrt(2)(|00〉 + |11〉)
}

operation PreparePhiMinus(register : Qubit[]) : Unit {
    H(register[0]);                 // |+0〉
    Z(register[0]);                 // |-0〉
    CNOT(register[0], register[1]); // 1/sqrt(2)(|00〉 - |11〉)
}

operation PreparePsiPlus(register : Qubit[]) : Unit {
    H(register[0]);                 // |+0〉
    X(register[1]);                 // |+1〉
    CNOT(register[0], register[1]); // 1/sqrt(2)(|01〉 + |10〉)
}


operation PreparePsiMinus(register : Qubit[]) : Unit {
    H(register[0]);                 // |+0〉
    Z(register[0]);                 // |-0〉
    X(register[1]);                 // |-1〉
    CNOT(register[0], register[1]); // 1/sqrt(2)(|01〉 - |10〉)
}
