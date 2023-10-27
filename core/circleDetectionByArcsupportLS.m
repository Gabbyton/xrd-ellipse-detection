function [circles, L,posi] = circleDetectionByArcsupportLS(I, Tac, Tni)
%input:
% I: input image
% Tac: angle of circular completeness
% Tni: ratio of support edge pixels on a circle
%output:
% circles: N by 3. (x,y,r)
% reference:
% 1��von Gioi R Grompone, Jeremie Jakubowicz, Jean-
% Michel Morel, and Gregory Randall, ��Lsd: a fast line
% segment detector with a false detection control.,�� IEEE
% transactions on pattern analysis and machine intelligence,
% vol. 32, no. 4, pp. 722�C732, 2010.
% 2��Truc Le and Ye Duan, ��Circle detection on images by
% line segment and circle completeness,�� in IEEE Inter-
% national Conference on Image Processing, 2016, pp.
% 3648�C3652.
    angleCoverage = Tac;%Բ�Ļ��ȸ���>=165��
    Tmin = Tni;%������ֵ0.5 
    unit_dis_tolerance = 2;%�ڵ��������̲�С��max(2,0.5%*minsize)
    normal_tolerance = 20;%�������̽Ƕ�20��
    t0 = clock;
    if(size(I,3)>1)
    [candidates,edge,normals] = generateCircleCandidates(I(:,:,3));%B����
    else
        [candidates,edge,normals] = generateCircleCandidates(I);
    end
%    figure;imshow(lsimg);
%    t1 = clock;
%    disp(['��ȡ��ѡԲ���ʱ��:',num2str(etime(t1,t0))]);
%    figure;imshow(edge);
    candidates = candidates';
    if(candidates(1) == 0)%��ʾû���ҵ���ѡԲ
        candidates =  zeros(0, 3);
    end
    posi = candidates;
    normals    = normals';
    [y, x]=find(edge);%�ҵ���0Ԫ�ص���(y)����(x)������
    [mylabels,labels, circles] = circleDetection(candidates,normals,[x, y], unit_dis_tolerance, normal_tolerance, Tmin, angleCoverage, edge);%���ĸ����� 0.5% 20�� 0.6 180�� 
    disp('-----------------------------------------------------------');
    disp(['running time:',num2str(etime(clock,t0)),'s']);
%     labels
%     size(labels)
%     size(y)
    warning('on', 'all');
     L = zeros(size(I, 1), size(I, 2));%����������ͼ��Iһ����С��0����L
     L(sub2ind(size(L), y, x)) = mylabels;%labels,���ȵ���edge_pixel_n x 1,�����i����Ե������ʶ���˵�j��Բ������б��Ϊj,����Ϊ0����С edge_pixel_n x 1;����ת���浽ͼ���У���ͼ���б��
%     figure;imshow(L==2);%LLL
%     imwrite((L==2),'D:\Graduate Design\��ͼ\edge_result.jpg');
end
%================================================================================================================================
%����1
%����
%points:     ��Ե���ص������(x,y),nx2,nΪ�ܹ��ı�Ե����
%lineLabels: ����Ӧ������(xi,yi)��ǣ���Ӧ������Ӧ���߶Σ�nx1,δ�����Ϊ0
%lines:      �߶β�����-B,A,xmid,ymid������(xmid,ymid)��Ӧ��Ӧ���߶��е㣬mx4��mΪ�ܹ�m���߶�
%���
%labels��    ���ȵ���n x 1,�����i����Ե������ʶ���˵�j��Բ������б��Ϊj,����Ϊ0����С n x 1
%circles��   ʶ�������Բ�ĺͰ뾶��ÿһ�и�ʽ��(x,y,r)
function [mylabels,labels, circles] = circleDetection(candidates, normals, points,  distance_tolerance, normal_tolerance, Tmin, angleCoverage,E)
    labels = zeros(size(points, 1), 1);
    mylabels = zeros(size(points, 1), 1);%����
    circles = zeros(0, 3);
    %�ͷź���Ϊ�յĲ���Ҳ��������candidates��СΪnCandidates x 3
    %candidates
    %angles = 0;%mycode 
    t1 = clock;
    angles = [340; 250; 160; 70];%�ǶȴӴ�С��֤��������
    angles(angles < angleCoverage) = [];%ֻ��������angleCoverage�Ĳ���
    if (isempty(angles) || angles(end) ~= angleCoverage)%���angelsΪ���ˣ�����angles��С��~=angleCoverage�����angleCoverage�������
        angles = [angles; angleCoverage];
    end
%    disp('��ʼ��һ��Բ�������ǶȽ�����֤����ʼangleLoop����ÿ��ѭ��������һ��angleCoverage���Ժ�ѡԲ������֤������Բ���ϵ��ڵ����ͨ�Է��������������������ȷ������Ӷ��ҵ���ЧԲ��ͬʱ�޳���ЧԲ');
    for angleLoop = 1 : length(angles)
        idx = find(labels == 0);%labels��СΪ��Ե��������enx1����ʼ��ʱlabelsȫΪ0���ҵ�labels�е���0������
        if (length(idx) < 2 * pi * (3 * distance_tolerance) * Tmin)%��idx����С��һ��ֵʱ
            break;
        end
        [L2,L, C, validCandidates] = subCircleDetection(points(idx, :), normals(idx, :), candidates, distance_tolerance, normal_tolerance, Tmin, angles(angleLoop),E,angleLoop);
        candidates = candidates(validCandidates, :);%����logical����validCandidates�����޳�����������Բ��ʣ�µ�Բ����������һ��angleloop��֤
      %  size(candidates)
      % disp(angleLoop)
        if (size(C, 1) > 0)
            for i = 1 : size(C, 1)
                flag = false;
                for j = 1 : size(circles, 1)
                    %��ʶ�������Բ���ܹ���֮ǰʶ�������Բ�ظ�
                    if (sqrt((C(i, 1) - circles(j, 1)) .^ 2 + (C(i, 2) - circles(j, 2)) .^ 2) <= distance_tolerance && abs(C(i, 3) - circles(j, 3)) <= distance_tolerance)
                        flag = true;
                        labels(idx(L == i)) = j;%����ظ��ˣ��ͰѸñ�ǩת�Ƶ�֮ǰ��Բ���档ע��ע�⣺idx�����������label�Ǳ�Ǹñ�Ե�������˵�j��Բ�ϣ�idx��labels����һά������(n x 1)��labels���Ե��points(n x 2)����һ��,��idx��һ��
                        %==================================================
                        mylabels(idx(L2 == i)) = j;
                        %==================================================
                        break;%������ѭ����������һ����ѭ��
                    end
                end
                if (~flag)%������ظ�������뵽ʶ���Բ(circles)��
                    labels(idx(L == i)) = size(circles, 1) + 1;
                    %=================================================================
                    %%��ʾ��ϳ�Բʱ���õ��ڵ�  my code 
                    mylabels(idx(L2 == i)) = size(circles, 1) + 1;%����
                    %=================================================================
                    circles = [circles; C(i, :)];
                end
            end
        end
    end
    t2 = clock;
%    disp(['�������֤ʱ�䣺',num2str(etime(t2,t1))]);
end
%================================================================================================================================
%����2
%����
%points:     ��Ե���ص������(x,y),nx2,nΪ�ܹ��ı�Ե����
%normals:    ÿһ����Ե���Ӧ���ݶ�������normals��СΪnx2����ʽΪ(xi,yi)
%list��      �����ѡ��Բ�ĺͰ뾶��ϣ�(x,y,r)����С candidate_n x 3. 
%���
%labels��    �����i����Ե�����ڼ�⵽�˵�j��Բ����labels��i�и�ֵΪj������Ϊ0.������pointsһ�£�n x 1
%circles:    �˴μ�⵽��Բ,(x,y,z),����⵽detectnum�������СΪdetectnum x 3
%validCandidates: list�ĺ�ѡԲ�У������i��Բ����⵽�˻��߲�����Բ����(Բ�����ڵ���������)�����i��λ��Ϊfalse(��ʼ��ʱΪtrue)����������һ��angleloop�ִ���֤ʱ�����޳���������Ҫ�ظ���֤��
%                 validCandidates�Ĵ�СΪ candidate_n x 1.
function [mylabels,labels, circles, validCandidates] = subCircleDetection(points, normals, list, distance_tolerance, normal_tolerance, Tmin, angleCoverage,E,angleLoop)
    labels = zeros(size(points, 1), 1);%��Ե���ص��������n,n x 1
    mylabels = zeros(size(points, 1), 1);%����
    circles = zeros(0, 3);
    maxRadius = min(max(points) - min(points));%���Ŀ��ܰ뾶(�˴��ɸ�Ϊ/2)
    validCandidates = true(size(list, 1), 1);%logical��������С candidate_n x 1
    convergence = list;
    for i = 1 : size(list, 1)
        circleCenter = list(i, 1 : 2);
        circleRadius = list(i, 3);
        tbins = min([180, floor(2 * pi * circleRadius * Tmin)]);%ѡ����
        circle_normals = points - repmat(circleCenter, size(points, 1), 1);
        circle_normals = circle_normals ./ repmat(sqrt(sum(circle_normals .^ 2, 2)), 1, 2);
        %�ҵ����Բlist(i,:)��Բ���ϵ��ڵ�,find�����жϵĽ��Ϊlogical���������ڵ����Ӧ��points�е��ж�Ӧλ��Ϊ1������Ϊ0����СΪn x 1
        %ͨ��find�ҳ���0Ԫ�ص�������inliers���Ǳ�������ӦΪ�ڵ���е�ֵ������ inlier_n x 1
        

        %���д��������⣬δ���м��Է����������ڽ���2 * distance_tolerance�ļ����෴�Ĵ����ڵ㱻���룬�Ӷ���ϴ���.
%           inliers = find(labels == 0 & abs(sqrt((points(:, 1) - circleCenter(1)) .^ 2 + (points(:, 2) - circleCenter(2)) .^ 2) - circleRadius) <= 2 * distance_tolerance & radtodeg(real(acos(abs(dot(normals, circle_normals, 2))))) <= normal_tolerance);
        %===============================================================================================================================================================
        %���д����Ϊ���´���,my code
         inliers = find(labels == 0 & abs(sqrt((points(:, 1) - circleCenter(1)) .^ 2 + (points(:, 2) - circleCenter(2)) .^ 2) - circleRadius) <= 2 * distance_tolerance);
         p_dot_temp = dot(normals(inliers,:), circle_normals(inliers,:), 2);
        p_cnt = sum(p_dot_temp>0);%����֮�٣���һ�μ���ͳ�ƣ�����ΪC����ʱע�����Բʱ�ڵ㼫�Ե�ѡȡ����
        if(p_cnt > size(inliers,1)*0.5)
            %��������,Ҳ�����ں����    
            inliers = inliers(p_dot_temp>0 & (radtodeg(real(acos(p_dot_temp))) <= normal_tolerance ));
        else
            %������ͬ,Ҳ�����ڰ���� 
            inliers = inliers(p_dot_temp<0 & (radtodeg(real(acos(-p_dot_temp))) <= normal_tolerance ));
        end
%           inliers = inliers((radtodeg(real(acos(abs(p_dot_temp)))) <= normal_tolerance ));
          
          inliers2 = inliers;
          inliers3 = 0;
        %=================================================================================================================================================================
        %��ͨ�������inliersΪ������ڱ�Ե������±�
%         size(points)
%          size(inliers)
%         size(points(inliers, :))
%         size(takeInliers(points(inliers, :), circleCenter, tbins))
        %��ͨ��������õ���Ч���ڵ㣬�ڵ��ᴿ��Ҳ����inliers�н�һ��������Ч��inliers,��������٣���Сinlier_n2 x 1��ע��ע�⣺inliers�д������points�е����±�
        inliers = inliers(takeInliers(points(inliers, :), circleCenter, tbins));
%         size(inliers)
        a = circleCenter(1); b = circleCenter(2); r = circleRadius; cnd = 0;
        [newa, newb, newr, newcnd] = fitCircle(points(inliers, :));
        
%         if angleLoop == 2   %mycode
%         dispimg = zeros(size(E,1),size(E,2),3);
%         dispimg(:,:,1) = E.*255;%��Ե��ȡ��������0-1ͼ��
%         dispimg(:,:,2) = E.*255;
%         dispimg(:,:,3) = E.*255;
%         for i = 1:length(inliers)
%         dispimg(points(inliers(i),2),points(inliers(i),1),:)=[0 0 255];
%         end
%         dispimg = drawCircle(dispimg,[newa newb],newr);
%         figure;
%         imshow(uint8(dispimg));
%         end

        if (newcnd == 0)%���������С���˷���ϵĶ��ó��Ľ��
            %��Բ�ĺ���Բ�ĵľ���С��4distance_tolerance����Ϊ����ϳ����Ĳ��ܺ�ԭ����Բ�Ĳ�ܶ�
            if (sqrt((newa - a) .^ 2 + (newb - b) .^ 2) <= 4 * distance_tolerance && abs(newr - r) <= 4 * distance_tolerance)
                circle_normals = points - repmat([newa, newb], size(points, 1), 1);
                circle_normals = circle_normals ./ repmat(sqrt(sum(circle_normals .^ 2, 2)), 1, 2);
                %������һ�����ڵ㣬��ͨ�Է������ڵ��ᴿ,��ε��µ��ڵ�����ں���������ȷ���
                newinliers = find(labels == 0 & abs(sqrt((points(:, 1) - newa) .^ 2 + (points(:, 2) - newb) .^ 2) - newr) <= 2 * distance_tolerance & radtodeg(real(acos(abs(dot(normals, circle_normals, 2))))) <= normal_tolerance);
                newinliers = newinliers(takeInliers(points(newinliers, :), [newa, newb], tbins));
                if (length(newinliers) >= length(inliers))
                    a = newa; b = newb; r = newr; cnd = newcnd;
                    inliers = newinliers;
                    inliers3 = newinliers;%my code��just test
                    %======================================================================
                    %�������
                    [newa, newb, newr, newcnd] = fitCircle(points(inliers, :));
                    if(newcnd == 0)
                        a = newa; b = newb; r = newr;
                    end
                    %=======================================================================
                end
            end
        end    
        
        if (length(inliers) >= floor(2 * pi * circleRadius * Tmin))%�ڵ���������Բ���ϵ�һ��������TminΪ������ֵ
            convergence(i, :) = [a, b, r];
            %��֮ǰ��Բ�ĺͰ뾶��������һ�£��ظ��ˣ���˰����Բ��̭(�����ͷ�ĺ����ظ���Բ��һ���ᱻ��̭)
            if (any(sqrt(sum((convergence(1 : i - 1, 1 : 2) - repmat([a, b], i - 1, 1)) .^ 2, 2)) <= distance_tolerance & abs(convergence(1 : i - 1, 3) - repmat(r, i - 1, 1)) <= distance_tolerance))
                validCandidates(i) = false;
            end
            %����ڵ���Բ��������angleCoverage��������
            completeOrNot =  isComplete(points(inliers, :), [a, b], tbins, angleCoverage);
            if (cnd == 0 && r < maxRadius && length(inliers) >= floor(2 * pi * r * Tmin) && completeOrNot)
                %�����������Բ��������distance_tolerance��Ҳ����ָ������Բ�ǲ�ͬ��
                if (all(sqrt(sum((circles(:, 1 : 2) - repmat([a, b], size(circles, 1), 1)) .^ 2, 2)) > distance_tolerance | abs(circles(:, 3) - repmat(r, size(circles, 1), 1)) > distance_tolerance))
                    %size(inliers)
                    line_normal = pca(points(inliers, :));%�õ�2x2��pca�任������˵ڶ��б������ڵ�ͳ�Ƴ����ݶ�
                    line_normal = line_normal(:, 2)';%ȡ���ڶ��в��ұ�Ϊ1 x 2 ��������
                    line_point = mean(points(inliers, :));%�ڵ�ȡƽ��
                    %��ֹ���ݵ���ڼ���
                    if (sum(abs(dot(points(inliers, :) - repmat(line_point, length(inliers), 1), repmat(line_normal, length(inliers), 1), 2)) <= distance_tolerance & radtodeg(real(acos(abs(dot(normals(inliers, :), repmat(line_normal, length(inliers), 1), 2))))) <= normal_tolerance) / length(inliers) < 0.8)
                         labels(inliers) = size(circles, 1) + 1;%��ǣ���Щ�ڵ��Ѿ��ù��ˣ��������¼�⵽Բ��
                         %==================================================================
                         if(all(inliers3) == 1)
                         mylabels(inliers3) = size(circles,1) + 1; %��ʾ��ϳ�Բʱ���õ��ڵ�  SSS
                         end
                         %==================================================================
                        circles = [circles; [a, b, r]];%����Բ���������ȥ
                        validCandidates(i) = false;%��i����ѡԲ������
                    end
                end
            end
        else
            validCandidates(i) = false;%�����������̭�ú�ѡԲ
        end
        
    end %for
end%fun
%================================================================================================================================
%����3
%�����߶ζ�L��Բ��
%���룺�߶ζ�AB,CD,L��С 1x8,�ֱ�洢��A(x1,y1),B(x2,y2),C(x3,y3),D(x4,y4)
%�����Բ��(x,y)
function center = initialEstimateCircleCenter(L)
    %center = ([L(1 : 2); L(5 : 6)] \ dot([L(3 : 4); L(7 : 8)], [L(1 : 2); L(5 : 6)], 2))';
    center = ([L(1 : 2); L(5 : 6)] \ dot([L(1 : 2); L(5 : 6)], [L(3 : 4); L(7 : 8)], 2))';
end
%����4
%Բ����С���˷����(�˴����Ը��ÿ���Բ��Ϸ���)
%���룺
%points: ��ͨ�Է�������ᴿ����ڵ㣬���СΪ fpn x 2,��ʽ(xi,yi)
%�����
%a   ����Ϻ��Բ�ĺ�����x
%b   ����Ϻ��Բ��������y
%c   ����Ϻ��Բ�İ뾶r
%cnd ��1��ʾ���ݴ��뷽�̺�������ģ�ֱ����ƽ��ֵ���ƣ�0��ʾ����������С���˷���ϵ�
function [a, b, r, cnd] = fitCircle(points)
%{
    A = [sum(points(:, 1)), sum(points(:, 2)), size(points, 1); sum(points(:, 1) .* points(:, 2)), sum(points(:, 2) .* points(:, 2)), sum(points(:, 2)); sum(points(:, 1) .* points(:, 1)), sum(points(:, 1) .* points(:, 2)), sum(points(:, 1))];
    %����С���˷�ʱ��A'A�����������ӽ�0������ζ�ŷ��������ԣ���ƽ��ֵ����
    if (abs(det(A)) < 1e-9)
        cnd = 1;
        a = mean(points(:, 1));
        b = mean(points(:, 2));
        r = min(max(points) - min(points));
        return;
    end
    cnd = 0;
    B = [-sum(points(:, 1) .* points(:, 1) + points(:, 2) .* points(:, 2)); -sum(points(:, 1) .* points(:, 1) .* points(:, 2) + points(:, 2) .* points(:, 2) .* points(:, 2)); -sum(points(:, 1) .* points(:, 1) .* points(:, 1) + points(:, 1) .* points(:, 2) .* points(:, 2))];
    t = A \ B;
    a = -0.5 * t(1);
    b = -0.5 * t(2);
    r = sqrt((t(1) .^ 2 + t(2) .^ 2) / 4 - t(3));
 %}
    A = [sum(points(:, 1) .* points(:, 1)),sum(points(:, 1) .* points(:, 2)),sum(points(:, 1)); sum(points(:, 1) .* points(:, 2)),sum(points(:, 2) .* points(:, 2)),sum(points(:, 2)); sum(points(:, 1)),sum(points(:, 2)),size(points, 1)]; 
    %����С���˷�ʱ��A'A�����������ӽ�0������ζ�ŷ��������ԣ���ƽ��ֵ����
    if (abs(det(A)) < 1e-9)
        cnd = 1;
        a = mean(points(:, 1));
        b = mean(points(:, 2));
        r = min(max(points) - min(points));
        return;
    end
    cnd = 0;
    B = [sum(-points(:, 1) .* points(:, 1) .* points(:, 1) - points(:, 1) .* points(:, 2) .* points(:, 2));sum(-points(:, 1) .* points(:, 1) .* points(:, 2) - points(:, 2) .* points(:, 2) .* points(:, 2)); sum(-points(:, 1) .* points(:, 1) - points(:, 2) .* points(:, 2))];
    t = A \ B;
    a = -0.5 * t(1);
    b = -0.5 * t(2);
    r = sqrt((t(1) .^ 2 + t(2) .^ 2) / 4 - t(3));
end
%================================================================================================================================
%����5
%����
%x     : ��ͨ�Է�������������2piRT���ᴿ����ڵ�(x,y)�������뵽�����ȷ�������.num x 2
%center: Բ��(x,y)  1 x 2
%tbins ����������
%angleCoverage: ��Ҫ�ﵽ��Բ������
%���
%result�� true or false����ʾ��Բ�����벻����
%longest_inliers:
function [result, longest_inliers] = isComplete(x, center, tbins, angleCoverage)
    [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%thetaΪ(-pi,pi)�ĽǶȣ�num x 1
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%theta�ĵ�i��Ԫ�����ڵ�j��bin����tt��i�б��Ϊj����Сnum x 1
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);
    longest_run = 0;
    start_idx = 1;
    end_idx = 1;
    while (start_idx <= tbins)
        if (h(start_idx) > 0)%�ҵ�bin��vote��һ������0��
            end_idx = start_idx;
            while (start_idx <= tbins && h(start_idx) > 0)%ֱ��bin��һ��С��0��
                start_idx = start_idx + 1;
            end
            inliers = [end_idx, start_idx - 1];%������Ϊ��ͨ����
            inliers = find(tt >= inliers(1) & tt <= inliers(2));%��tt���ҵ����ڴ�������ڵ������
            run = max(theta(inliers)) - min(theta(inliers));%�ǶȲ�
            if (longest_run < run)%�˾���Ϊ���ҵ���������������ͨ�Ŀ��
                longest_run = run;
                longest_inliers = inliers;
            end
        end
        start_idx = start_idx + 1;
    end
    if (h(1) > 0 && h(tbins) > 0)%�����һ��bin�����һ��bin������0���п��������ͨ������ͷβ�������������
        start_idx = 1;
        while (start_idx < tbins && h(start_idx) > 0)%�ҵ�bin��vote��һ������0��
            start_idx = start_idx + 1;
        end
        end_idx = tbins;%end_idxֱ�Ӵ���β����ʼ������
        while (end_idx > 1 && end_idx > start_idx && h(end_idx) > 0)
            end_idx = end_idx - 1;
        end
        inliers = [start_idx - 1, end_idx + 1];
        run = max(theta(tt <= inliers(1)) + 2 * pi) - min(theta(tt >= inliers(2)));
        inliers = find(tt <= inliers(1) | tt >= inliers(2));
        if (longest_run < run)
            longest_run = run;
            longest_inliers = inliers;
        end
    end
    %������ͨ�Ŀ�ȴ�����angleCoverage��������Ȼ�����ͨ���С�ڣ����������㹻��
    longest_run_deg = radtodeg(longest_run);
    h_greatthanzero_num = sum(h>0);
    result =  longest_run_deg >= angleCoverage || h_greatthanzero_num * (360 / tbins) >= min([360, 1.2 * angleCoverage]);
end
%================================================================================================================================
%����6
%��ͨ�Է�������Բ���ϵ��ڵ�����ᴿ
%����
%x��Բ���ϵ��ڵ�(x,y),��Ϊinlier_n x 2 
%center��һ��Բ��Բ��(x,y) 1x2
%tbins: ���� = min(180,2pi*R*T) 
%���
%idx��Ϊ��xһ�����ģ�inlier_n x 1��logical������������Ч������һ����ͨ���ȵ��ڵ㣬��Ӧλ����Ч��Ϊ1������Ϊ0
function idx = takeInliers(x, center, tbins)
   [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%�õ�[-pi,pi]�ķ�λ�ǣ��ȼ��� theta = atan2(x(:, 2) - center(2) , x(:, 1) - center(1)); 
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%���ڵ������[1 tbins]
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);%hΪֱ��ͼ[1 tbins]��ͳ�ƽ��
    mark = zeros(tbins, 1);
    compSize = zeros(tbins, 1);
    nComps = 0;
    queue = zeros(tbins, 1);
    du = [-1, 1];
    for i = 1 : tbins
        if (h(i) > 0 && mark(i) == 0)%������ڵ�i�������ڵ�ֵ����0����mark(i)Ϊ0
            nComps = nComps + 1;
            mark(i) = nComps;%��ǵ�nComps����ͨ����
            front = 1; rear = 1;
            queue(front) = i;%���÷���������У����Դ˿�ʼ����
            while (front <= rear)
                u = queue(front);
                front = front + 1;
                for j = 1 : 2
                    v = u + du(j);
                    if (v == 0)
                        v = tbins;
                    end
                    if (v > tbins)
                        v = 1;
                    end
                    if (mark(v) == 0 && h(v) > 0)
                        rear = rear + 1;
                        queue(rear) = v;
                        mark(v) = nComps;%��ǵ�nComps����ͨ����
                    end
                end
            end
            compSize(nComps) = sum(ismember(tt, find(mark == nComps)));%�õ�������ͨ��ΪnComps���ڵ�����
        end
    end
    compSize(nComps + 1 : end) = [];
    maxCompSize = max(compSize);
    validComps = find(compSize >= maxCompSize * 0.1 & compSize > 10);%���ڵ��������ͨ���ȵ�0.1������ͨ��������Ч��
    validBins = find(ismember(mark, validComps));%��Ч�ķ���
    idx = ismember(tt, validBins);%��Ч���ڵ�
end
