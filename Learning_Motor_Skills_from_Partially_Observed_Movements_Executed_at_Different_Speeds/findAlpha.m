% find alpha for given trajectory, such that z = alpha*t and 0 <= z <= 100

function alpha = findAlpha(trajectory)

    length_trajectory = size(trajectory,1);
    alpha = 100/length_trajectory;

end