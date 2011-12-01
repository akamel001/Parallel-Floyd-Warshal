function D = floyd3(A)

  %D = floydinit(A);
  n = length(A);
  b = 5;

  for kk = 1:b:n
    kf = min(n,kk+b-1);

    % Update diagonal block for b steps
    for k = kk:kf
      for j = 1:kf
        for i = 1:kf
          D(i,j) = min(D(i,j), D(i,k) + D(k,j));
        end
      end
    end

    % Update leading columns for b steps
    for k = kk:kf
      for j = 1:kf
        for i = kf+1:n
          D(i,j) = min(D(i,j), D(i,k) + D(k,j));
        end
      end
    end

    % Update trailing matrix for b steps
    for k = kk:kf
      for j = kf+1:n
        for i = 1:n
          D(i,j) = min(D(i,j), D(i,k) + D(k,j));
        end
      end
    end

  end
