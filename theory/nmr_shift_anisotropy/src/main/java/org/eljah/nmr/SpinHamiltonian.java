package org.eljah.nmr;

import java.util.Arrays;

/** Exact Hermitian 4-by-4 I=3/2 Hamiltonian; all frequency quantities are MHz. */
final class SpinHamiltonian {
    private static final double ROOT3 = Math.sqrt(3);
    private static final double[] M = {1.5, .5, -.5, -1.5};
    private static final double[] LADDER_HALF = {ROOT3 / 2, 1, ROOT3 / 2};
    final double[] eigenvalues = new double[4];
    final int[] order = {0, 1, 2, 3};
    private final double[] ar = new double[16], ai = new double[16];
    private final double[] vr = new double[16], vi = new double[16];
    private final double[] tx = new double[2], ty = new double[2], tz = new double[2];

    void solve(double zx, double zy, double zz, double nuQ, double eta) {
        Arrays.fill(ar, 0); Arrays.fill(ai, 0);
        Arrays.fill(vr, 0); Arrays.fill(vi, 0);
        for (int i = 0; i < 4; i++) {
            ar[5 * i] = -zz * M[i] + (i == 0 || i == 3 ? nuQ / 2 : -nuQ / 2);
            vr[5 * i] = 1;
        }
        for (int i = 0; i < 3; i++) {
            int a = 4 * i + i + 1, b = 4 * (i + 1) + i;
            ar[a] = ar[b] = -zx * LADDER_HALF[i];
            ai[a] = zy * LADDER_HALF[i]; ai[b] = -ai[a];
        }
        ar[2] = ar[8] = ar[7] = ar[13] = nuQ * eta * ROOT3 / 6;
        diagonalize();
    }

    private void diagonalize() {
        for (int sweep = 0; sweep < 30; sweep++) {
            double largest = 0;
            for (int p = 0; p < 3; p++) for (int q = p + 1; q < 4; q++) {
                int pq = 4 * p + q;
                double magnitude = Math.hypot(ar[pq], ai[pq]);
                largest = Math.max(largest, magnitude);
                if (magnitude < 1e-13) continue;
                double ur = ar[pq] / magnitude, ui = ai[pq] / magnitude;
                double tau = (ar[5 * q] - ar[5 * p]) / (2 * magnitude);
                double t = Math.copySign(1.0, tau) / (Math.abs(tau) + Math.hypot(1, tau));
                double c = 1 / Math.sqrt(1 + t * t), s = t * c;
                ar[5 * p] -= t * magnitude;
                ar[5 * q] += t * magnitude;
                ar[pq] = ar[4 * q + p] = 0;
                ai[pq] = ai[4 * q + p] = 0;
                for (int k = 0; k < 4; k++) {
                    if (k != p && k != q) {
                        int kp = 4 * k + p, kq = 4 * k + q;
                        double pr = ar[kp], pi = ai[kp], qr = ar[kq], qi = ai[kq];
                        ar[kp] = c * pr - s * (ur * qr + ui * qi);
                        ai[kp] = c * pi - s * (ur * qi - ui * qr);
                        ar[kq] = s * (ur * pr - ui * pi) + c * qr;
                        ai[kq] = s * (ur * pi + ui * pr) + c * qi;
                        ar[4 * p + k] = ar[kp]; ai[4 * p + k] = -ai[kp];
                        ar[4 * q + k] = ar[kq]; ai[4 * q + k] = -ai[kq];
                    }
                    int kp = 4 * k + p, kq = 4 * k + q;
                    double pr = vr[kp], pi = vi[kp], qr = vr[kq], qi = vi[kq];
                    vr[kp] = c * pr - s * (ur * qr + ui * qi);
                    vi[kp] = c * pi - s * (ur * qi - ui * qr);
                    vr[kq] = s * (ur * pr - ui * pi) + c * qr;
                    vi[kq] = s * (ur * pi + ui * pr) + c * qi;
                }
            }
            if (largest < 1e-12) break;
            if (sweep == 29) throw new ArithmeticException("Hermitian Jacobi did not converge.");
        }
        for (int i = 0; i < 4; i++) order[i] = i;
        for (int i = 1; i < 4; i++) {
            int v = order[i], j = i;
            while (j > 0 && ar[5 * order[j - 1]] > ar[5 * v]) { order[j] = order[j - 1]; j--; }
            order[j] = v;
        }
        for (int i = 0; i < 4; i++) eigenvalues[i] = ar[5 * order[i]];
    }

    /** Archive RF weights, or optional full axial G.I RF weights. */
    double weight(int levelA, int levelB, double nx, double ny, double nz,
                  double gPerpendicular, double gParallel) {
        int a = order[levelA], b = order[levelB];
        double xr = 0, xi = 0, yr = 0, yi = 0, zr = 0, zi = 0;
        for (int k = 0; k < 4; k++) {
            int ka = 4 * k + a, kb = 4 * k + b;
            zr += M[k] * (vr[ka] * vr[kb] + vi[ka] * vi[kb]);
            zi += M[k] * (vr[ka] * vi[kb] - vi[ka] * vr[kb]);
        }
        for (int k = 0; k < 3; k++) {
            int a0 = 4 * k + a, a1 = 4 * (k + 1) + a;
            int b0 = 4 * k + b, b1 = 4 * (k + 1) + b;
            double pr = vr[a0] * vr[b1] + vi[a0] * vi[b1];
            double pi = vr[a0] * vi[b1] - vi[a0] * vr[b1];
            double qr = vr[a1] * vr[b0] + vi[a1] * vi[b0];
            double qi = vr[a1] * vi[b0] - vi[a1] * vr[b0];
            double l = LADDER_HALF[k];
            xr += l * (pr + qr); xi += l * (pi + qi);
            yr += l * (pi - qi); yi += l * (qr - pr);
        }
        xr *= gPerpendicular; xi *= gPerpendicular;
        yr *= gPerpendicular; yi *= gPerpendicular;
        zr *= gParallel; zi *= gParallel;
        double parallelR = nx * xr + ny * yr + nz * zr;
        double parallelI = nx * xi + ny * yi + nz * zi;
        double result = (xr * xr + xi * xi + yr * yr + yi * yi + zr * zr + zi * zi
            - parallelR * parallelR - parallelI * parallelI) / 2;
        return Math.max(0, result);
    }

    double orthogonalityError() {
        double largest = 0;
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
            double r = 0, im = 0;
            for (int k = 0; k < 4; k++) {
                int ka = 4 * k + a, kb = 4 * k + b;
                r += vr[ka] * vr[kb] + vi[ka] * vi[kb];
                im += vr[ka] * vi[kb] - vi[ka] * vr[kb];
            }
            largest = Math.max(largest, Math.hypot(r - (a == b ? 1 : 0), im));
        }
        return largest;
    }
}
