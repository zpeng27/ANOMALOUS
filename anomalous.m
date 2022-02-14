function [R] = anomalous(X, A, L, alpha, beta, gamma, phi,niters)

    [d,n] = size(X);
%     Dr = diag(1+rand(1,n));
%     Dw1 = diag(1+rand(1,n));
%     Dw2 = diag(1+rand(1,d));
    Dr = eye(n);
    Dw1 = eye(n);
    Dw2 = eye(d);
%     Dr = diag(ones(n,1)*1.5);
%     Dw1 = diag(ones(n,1)*1.5);
%     Dw2 = diag(ones(d,1)*1.5);
    R = X*inv(eye(n)+gamma*Dr+phi*L);
    [P,theta1] = eig(X'*X,'nobalance');
    [Q,theta2] = eig(X*X','nobalance');
    P = gs(P);
    Q = gs(Q);
    
    for iter = 1:niters
        %% update W
        H = X'*X*X'-X'*R*X';
        Dw_star = zeros(n,d);
        for i = 1:d
            Dw_star((1:n),i) = Dw2(i,i);
        end
        temp1 = P'*H*Q;
        M = P'*Dw1*P;
        Y = zeros(n,d);
        for i = 1:d
            temp2 = theta1.*diag(ones(n,1)*theta2(i,i))+alpha*M+beta*eye(n).*diag(Dw_star(:,i));
            Y(:,i) = inv(temp2)*temp1(:,i);
        end
        W = P*Y*Q';
        Dw1 = diag(0.5./sqrt(sum((W.^2),2)+eps));
        Dw2 = diag(0.5./sqrt(sum(((W').^2),2)+eps));
        
        %% update R
        R = (X-X*W*X)*inv(eye(n)+gamma*Dr+phi*L);
        Dr = diag(0.5./sqrt(sum(((R').^2),2)+eps));
        
        %% check if the objective function converges
        obj(iter) = norm(X-X*W*X-R,'fro')^2 + alpha*sum(sqrt(sum((W.^2),2))) + beta*sum(sqrt(sum(((W').^2),2))) + gamma*sum(sqrt(sum(((R').^2),2))) + phi*trace(R*L*R');
        fprintf('the object value in iter %d is %f\n', iter, obj(iter));
        if (iter >= 2 && abs(obj(iter)-obj(iter-1))<1e-3)
            break;
        end
    end
end