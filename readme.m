function [cost] = safeCostFunctionADeltav2(u,f,status0,Np,Nc,punish_weights,risk_th)
%safeCostFunction 非线性代价函数，包括与等高线的贴合程度、动力学操纵不突变、末态安全
% 由决策层做出需要做出换道的指令，这里先默认按照向右换道写，所以令末态的y=3.5*0.5,delta=0
% Np = Nc现在控制域Nc=2,预测域相同Np=25
% Nc = 2;

dt = 0.1;
T = dt;
l=2.5;
% 复杂度太高，不太可行
% shape: 3T*1 = 3T*3 * 3*1 + 3T*2T * 2T*1
% status_predict = A*status0+B*u;

X=status0(1,1);
Y=status0(2,1);
V=status0(3,1);
PHI=status0(4,1);


X_predict = zeros(Np,1);
Y_predict = zeros(Np,1);
V_predict = zeros(Np,1);
PHI_predict = zeros(Np,1);

a = zeros(Np,1); %控制量a
delta_f = zeros(Np,1); %控制量delta_f

%状态更新
for i =1:1:Np
    if i == 1
        a(i) = u(1);
        delta_f(i,1) = u(2);
        V_predict(i) = V + a(i,1)*T;
        X_predict(i) = X +T*V_predict(i)*cos(PHI);
        Y_predict(i) = Y +T*V_predict(i)*sin(PHI);
        PHI_predict(i) = PHI+T*V_predict(i)*tan(delta_f(i))/l;
    else
        if i<=Nc
            a(i,1)=u(2*i-1);
            delta_f(i,1)=u(2*i);
        else
            a(i,1)=u(2*Nc-1);
            delta_f(i,1)=0;
        end
        V_predict(i) = V_predict(i-1) + a(i)*T;

        X_predict(i)=X_predict(i-1)+T*V_predict(i)*cos(PHI_predict(i-1));
        Y_predict(i)=Y_predict(i-1)+T*V_predict(i)*sin(PHI_predict(i-1));
        PHI_predict(i)=PHI_predict(i-1)+T*V_predict(i)*tan(delta_f(i))/l;
    end

end



a_ctrl = u(1:2:2*Nc);
delta_ctrl = u(2:2:2*Nc);

a_ctrl = a_ctrl';
delta_ctrl = delta_ctrl';

% 行驶过轨迹的安全场的值,f_map precision is 0.05;
f_ix = int64(X_predict/0.05)+1;
f_iy = int64(Y_predict/0.05)+1;
% 简单保证不会超出地图
f_iy(f_iy>166) = 166;
f_iy(f_iy<1)=1;f_ix(f_ix<1)=1;
f_ix(f_ix>4000) = 4000;

f_iTime = 1:Np;

f_idx = sub2ind(size(f), f_iy, f_ix, f_iTime');


f_predict = f(f_idx);


qf = punish_weights.f;
qa = eye(Nc)*punish_weights.a;
qdelta = eye(Nc)*punish_weights.delta;
qy = eye(Np)*punish_weights.y;
qPhi = punish_weights.phiEnd;


costf = qf*(f_predict(:)-risk_th)'*(f_predict(:)-risk_th);
costf2 = qf*(f_predict(:))'*(f_predict(:));
costa = a_ctrl'*qa*a_ctrl;
costDelta = delta_ctrl'*qdelta*delta_ctrl;
% costY = (Y_predict - 1.75)'*qy*(Y_predict - 1.75);
costY = (Y_predict - 5.25)'*qy*(Y_predict - 5.25);


costPhi = qPhi*PHI_predict(end,1)^2;

% costs = [costf costa costDelta costY costPhi];
cost = costf+costa+costDelta+costY+costPhi+costf2;

end

