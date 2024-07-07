function result = gcd_of_three(a, b, c, threshold)
    min_value = min([a, b, c]);
    result = min_value;
    for i = threshold:-1:1
        if mod(a, i) == 0 && mod(b, i) == 0 && mod(c, i) == 0
            result = i;
            break;  % Exit loop once we find the GCD
        end
    end
end
