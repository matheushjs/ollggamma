context("Odd Log-Logistic Generalized Gamma distribution");

require(ollggamma);

test_that("dollggamma integrates to 1", {
	result = integrate(function(x) dollggamma(x, 1, 1, 1, 1), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dollggamma(x, 0.8, 0.2, 1, 0.5), subdivisions=1000000, 0, 1000000)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dollggamma(x, 0.5, 0.5, 0.2, 1.2), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dollggamma(x, 0.5, 5, 2, 1), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dollggamma(x, 0.5, 1, 2, 1), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);

	result = integrate(function(x) dollggamma(x, 5, 1, 0.5, 1), 0, Inf)$value;
	expect_equal(result, 1, tolerance=0.001);
});

test_that("pollggamma correctly reflects dollggamma", {
	expected = integrate(function(x) dollggamma(x, 1, 1, 1, 1), 0, 10)$value;
	result   = pollggamma(10, 1, 1, 1, 1);
	expect_equal(result, expected);

	result   = pollggamma(1e-10, 1, 1, 1, 1.5);
	expect_equal(result, 0);

	result   = pollggamma(20, 1, 1, 1, 1.5);
	expect_equal(result, 1);

	expected = integrate(function(x) dollggamma(x, 0.8, 0.5, 1, 1.2), 0, 10)$value;
	result   = pollggamma(10, 0.8, 0.5, 1, 1.2);
	expect_equal(result, expected, tolerance=1e-5);

	result   = pollggamma(0, 0.8, 0.5, 1, 1.2);
	expect_equal(result, 0);

	result   = pollggamma(1000, 0.8, 0.5, 1, 1.2);
	expect_equal(result, 1);

	expected = integrate(function(x) dollggamma(x, 1.2, 0.5, 1.5, 0.5), 0, 10)$value;
	result   = pollggamma(10, 1.2, 0.5, 1.5, 0.5);
	expect_equal(result, expected, tolerance=1e-5);

	result   = pollggamma(0, 1.2, 0.5, 1.5, 0.5);
	expect_equal(result, 0);

	result   = pollggamma(10000, 1.2, 0.5, 1.5, 0.5);
	expect_equal(result, 1);
});

test_that("qollggamma correctly inverts pollggamma", {
	a = 1; b = 1; k = 1; lambda = 1;
	expected = c(1, 10);
	expect_equal(qollggamma(pollggamma(expected, a, b, k, lambda), a, b, k, lambda), expected);

	a = 0.5; b = 1.5; k = 1.5; lambda = 1.2;
	expected = c(0.5, 1, 3);
	expect_equal(qollggamma(pollggamma(expected, a, b, k, lambda), a, b, k, lambda), expected);

	a = 1.5; b = 0.5; k = 3; lambda = 0.5;
	expected = c(0.5, 1, 3);
	expect_equal(qollggamma(pollggamma(expected, a, b, k, lambda), a, b, k, lambda), expected);
});

test_that("random number generation", {
	set.seed(72);

	a = 1; b = 1; k = 1; lambda = 1;
	r = rollggamma(100000, a, b, k, lambda);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pollggamma(x, a, b, k, lambda)));
	expect_lte(maxError, 1e-2);

	a = 0.3; b = 1.5; k = 1.5; lambda = 1.2;
	r = rollggamma(100000, a, b, k, lambda);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pollggamma(x, a, b, k, lambda)));
	expect_lte(maxError, 1e-2);

	a = 1.5; b = 5; k = 0.3; lambda = 0.5;
	r = rollggamma(100000, a, b, k, lambda);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pollggamma(x, a, b, k, lambda)));
	expect_lte(maxError, 1e-2);

	a = 3; b = 1.8; k = 0.5; lambda = 1.5;
	r = rollggamma(100000, a, b, k, lambda);
	e = ecdf(r);
	x = seq(0.001, 10, length=10000);
	maxError = max(abs(e(x) - pollggamma(x, a, b, k, lambda)));
	expect_lte(maxError, 1e-2);

	set.seed(NULL);
});
