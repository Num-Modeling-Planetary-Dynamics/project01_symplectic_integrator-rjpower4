kepler_error(ecc, ma, ea) = ea - ecc * sin(ea) - ma 
kepler_error_d1(ecc, ma, ea) = 1 - ecc * cos(ea)
kepler_error_d2(ecc, ma, ea) = ecc * sin(ea)
kepler_error_d3(ecc, ma, ea) = ecc * cos(ea)

function solve_kepler(updater, ecc, ma; tolerance=1e-12, max_iter=20)
    guess = ma
    error = kepler_error(ecc, ma, guess)

    iter = 1
    while abs(error) > tolerance && iter <= max_iter
        guess = updater(ecc, ma, guess)
        error = kepler_error(ecc, ma, guess)
        iter += 1
    end

    
    return guess
end