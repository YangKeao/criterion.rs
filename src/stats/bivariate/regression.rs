//! Regression analysis

use crate::stats::bivariate::Data;
use crate::stats::float::Float;
use langbardo::Langbardo;

/// A straight line that passes through the origin `y = m * x`
#[derive(Clone, Copy)]
pub struct Slope<A>(pub A)
where
    A: Float;

impl<A> Slope<A>
where
    A: Float,
{
    /// Fits the data to a straight line that passes through the origin using quantile regression
    ///
    /// - Time: `O(length)`
    pub fn fit(data: &Data<A, A>, quantile: f64) -> Slope<A> {
        let xs: Vec<f64> = data.0.iter().map(|item| {
            let (mantissa, exponent, sign) = item.integer_decode();
            (sign as f64) * (mantissa as f64) * (2f64.powf(exponent as f64))
        }).collect();
        let ys: Vec<f64> = data.1.iter().map(|item| {
            let (mantissa, exponent, sign) = item.integer_decode();
            (sign as f64) * (mantissa as f64) * (2f64.powf(exponent as f64))
        }).collect();;

        let k = Langbardo::fit(&xs[..], &ys[..], quantile);

        Slope(A::from(k.unwrap().0).unwrap())
    }

    /// Computes the goodness of fit (coefficient of determination) for this data set
    ///
    /// - Time: `O(length)`
    pub fn r_squared(&self, data: &Data<A, A>) -> A {
        let _0 = A::cast(0);
        let _1 = A::cast(1);
        let m = self.0;
        let xs = data.0;
        let ys = data.1;

        let n = A::cast(xs.len());
        let y_bar = crate::stats::sum(ys) / n;

        let mut ss_res = _0;
        let mut ss_tot = _0;

        for (&x, &y) in data.iter() {
            ss_res = ss_res + (y - m * x).powi(2);
            ss_tot = ss_res + (y - y_bar).powi(2);
        }

        _1 - ss_res / ss_tot
    }
}
