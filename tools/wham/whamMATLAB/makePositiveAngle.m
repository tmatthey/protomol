function result = makePositiveAngle(angle)
    % Convert a vector of angles in the range [-PI PI] to
    % the range of [0 2*PI]
    adjustment = 2 * pi * (angle < 0);
    result = adjustment + angle;