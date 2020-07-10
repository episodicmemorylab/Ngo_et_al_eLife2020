function [X_y, Y_x, xy_corr, yx_corr]= mt_pwrcorr_ortho(method,X,Y,rotate)
% Orthogonalization of complex (time-frequency) time series using the
% method from Hipp et all, Nat Neurosci 2011. 
%
%Input:
% method  - 'orthoregress' simple orthogonalization using global regression
%           'hipp' Hipp et all, Nat Neurosci 2011, pointwise
%           orthogonalization
% X,Y     - complex input signals that are to be orthogonalized wrt each
%           other
% rotate  - for Hipp method: correctly rotate the complex coefficients (not 
%           necessary if you
%           only want to get envelope correlation)
%
%Returns:
% X_y     - X orthogonalized wrt Y
% Y_x     - Y orthogonalized wrt X
% xy_corr - envelope correlation between X_y and Y
% yx_corr - envelope correlation between Y_x and X

if nargin<4, rotate=0; end


switch(method)
    case 'hipp'
        % Orientation in complex plane orthogonal to x / y (see Hipp
        % supplemental material)

        if rotate
            % !! This just rotates the phase and does not need to be done 
            % here if we just want the envelope correlation !!
            e_x= 1i*X ./ abs(X);
            e_y= 1i*Y ./ abs(Y);

            Y_x= imag(Y .* conj(X) ./ abs(X)) .* e_x;
            X_y= imag(X .* conj(Y) ./ abs(Y)) .* e_y;
        else
            % Orthonalize
            Y_x= imag(Y .* conj(X) ./ abs(X));
            X_y= imag(X .* conj(Y) ./ abs(Y));
        end
    case 'orthoregress'
        B= regress(real(X)', real(Y)');
        X_y = X- B* Y;

        B= regress(real(Y)', real(X)');
        Y_x = Y- B* X;
    otherwise
        error('Unknown method %s', method)
end

if nargout>2
    % Determine envelope correlations
    xy_corr= corr(abs(X)', abs(Y_x)');
    yx_corr= corr(abs(Y)', abs(X_y)');
end
 