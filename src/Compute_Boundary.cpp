//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#include <cmath>

#include "Blunt_Solver.h"
#include "Geometry.h"
#include "QlQr_Solver.h"

void Compute_Boundary_Blunt_Body() {
    VInt2D& marker = mesh->Get_Marker_Q();
    int ist, ied, jst, jed;
    Get_IJK_Region(ist, ied, jst, jed);

    // 上边界：outflow
    for (int i = ist; i <= ied; i++) {
        for (int j = jed; j <= jed + 2; j++) {
            qField[i][j][IR] = qField[i][j - 1][IR];
            qField[i][j][IU] = qField[i][j - 1][IU];
            qField[i][j][IV] = qField[i][j - 1][IV];
            qField[i][j][IP] = qField[i][j - 1][IP];
        }
    }

    // 下边界：outflow
    for (int i = ist; i <= ied; i++) {
        for (int j = jst; j >= 0; j--) {
            qField[i][j][IR] = qField[i][j + 1][IR];
            qField[i][j][IU] = qField[i][j + 1][IU];
            qField[i][j][IV] = qField[i][j + 1][IV];
            qField[i][j][IP] = qField[i][j + 1][IP];
        }
    }

    // 右边界：outflow
    for (int j = jst; j <= jed; j++) {
        for (int i = ied; i <= ied + 2; i++) {
            if (marker[i][j] == 0)
                continue;

            qField[i][j][IR] = qField[i - 1][j][IR];
            qField[i][j][IU] = qField[i - 1][j][IU];
            qField[i][j][IV] = qField[i - 1][j][IV];
            qField[i][j][IP] = qField[i - 1][j][IP];
        }
    }
    // 左边界：超声速入口
    for (int i = 0; i <= ist; i++) {
        for (int j = jst; j <= jed; j++) {
            qField[i][j][IR] = 1.0;
            qField[i][j][IU] = 3.0;
            qField[i][j][IV] = 0.0;
            qField[i][j][IP] = 0.71429;
        }
    }

    // 上物面, slip wall
    int JJ = jst + Jw2;
    for (int i = ist + Iw; i <= ied; i++) {
        qField[i][JJ][IV] = 0.0;

        qField[i][JJ - 1][IR] = qField[i][JJ + 1][IR];
        qField[i][JJ - 1][IU] = qField[i][JJ + 1][IU];
        qField[i][JJ - 1][IV] = -qField[i][JJ + 1][IV];
        qField[i][JJ - 1][IP] = qField[i][JJ + 1][IP];

        qField[i][JJ - 2][IR] = qField[i][JJ + 2][IR];
        qField[i][JJ - 2][IU] = qField[i][JJ + 2][IU];
        qField[i][JJ - 2][IV] = -qField[i][JJ + 2][IV];
        qField[i][JJ - 2][IP] = qField[i][JJ + 2][IP];
    }

    // 下物面, slip wall
    JJ = jst + Jw1;
    for (int i = ist + Iw; i <= ied; i++) {
        qField[i][JJ][IV] = 0.0;

        qField[i][JJ + 1][IR] = qField[i][JJ - 1][IR];
        qField[i][JJ + 1][IU] = qField[i][JJ - 1][IU];
        qField[i][JJ + 1][IV] = -qField[i][JJ - 1][IV];
        qField[i][JJ + 1][IP] = qField[i][JJ - 1][IP];

        qField[i][JJ + 2][IR] = qField[i][JJ - 2][IR];
        qField[i][JJ + 2][IU] = qField[i][JJ - 2][IU];
        qField[i][JJ + 2][IV] = -qField[i][JJ - 2][IV];
        qField[i][JJ + 2][IP] = qField[i][JJ - 2][IP];
    }

    // 左物面, slip wall
    int II = ist + Iw;
    for (int j = jst + Jw1; j <= jst + Jw2; j++) {
        qField[II][j][IU] = 0.0;

        qField[II + 1][j][IR] = qField[II - 1][j][IR];
        qField[II + 1][j][IU] = -qField[II - 1][j][IU];
        qField[II + 1][j][IV] = qField[II - 1][j][IV];
        qField[II + 1][j][IP] = qField[II - 1][j][IP];

        qField[II + 2][j][IR] = qField[II - 2][j][IR];
        qField[II + 2][j][IU] = -qField[II - 2][j][IU];
        qField[II + 2][j][IV] = qField[II - 2][j][IV];
        qField[II + 2][j][IP] = qField[II - 2][j][IP];
    }

    // 角点处理，左上，设置参考角点标号
    II = ist + Iw;
    JJ = jst + Jw2;
    // 1点，对角方向
    qField[II + 2][JJ - 2][IR] = qField[II][JJ][IR];
    qField[II + 2][JJ - 2][IU] = -qField[II][JJ][IU];
    qField[II + 2][JJ - 2][IV] = -qField[II][JJ][IV];
    qField[II + 2][JJ - 2][IP] = qField[II][JJ][IP];

    // 2点, x方向取值
    qField[II + 1][JJ - 2][IR] = qField[II - 1][JJ - 2][IR];
    qField[II + 1][JJ - 2][IU] = -qField[II - 1][JJ - 2][IU];
    qField[II + 1][JJ - 2][IV] = qField[II - 1][JJ - 2][IV];
    qField[II + 1][JJ - 2][IP] = qField[II - 1][JJ - 2][IP];

    // 3点, y方向取值
    qField[II + 2][JJ - 1][IR] = qField[II + 2][JJ + 1][IR];
    qField[II + 2][JJ - 1][IU] = qField[II + 2][JJ + 1][IU];
    qField[II + 2][JJ - 1][IV] = -qField[II + 2][JJ + 1][IV];
    qField[II + 2][JJ - 1][IP] = qField[II + 2][JJ + 1][IP];

    // 4点， x方向取值和y方向取值的算术平均
    qField[II + 1][JJ - 1][IR] = 0.5 * (qField[II - 1][JJ - 1][IR] + qField[II + 1][JJ + 1][IR]);
    qField[II + 1][JJ - 1][IU] = 0.5 * (qField[II - 1][JJ - 1][IU] + qField[II + 1][JJ + 1][IU]);
    qField[II + 1][JJ - 1][IV] = 0.5 * (qField[II - 1][JJ - 1][IV] + qField[II + 1][JJ + 1][IV]);
    qField[II + 1][JJ - 1][IP] = 0.5 * (qField[II - 1][JJ - 1][IP] + qField[II + 1][JJ + 1][IP]);

    //=========================
    // 角点处理，左下，设置参考角点标号
    II = ist + Iw;
    JJ = jst + Jw1;

    // 1点，对角方向
    qField[II + 2][JJ + 2][IR] = qField[II][JJ][IR];
    qField[II + 2][JJ + 2][IU] = -qField[II][JJ][IU];
    qField[II + 2][JJ + 2][IV] = -qField[II][JJ][IV];
    qField[II + 2][JJ + 2][IP] = qField[II][JJ][IP];

    // 2点, x方向取值
    qField[II + 1][JJ + 2][IR] = qField[II - 1][JJ + 2][IR];
    qField[II + 1][JJ + 2][IU] = -qField[II - 1][JJ + 2][IU];
    qField[II + 1][JJ + 2][IV] = qField[II - 1][JJ + 2][IV];
    qField[II + 1][JJ + 2][IP] = qField[II - 1][JJ + 2][IP];

    // 3点, y方向取值
    qField[II + 2][JJ + 1][IR] = qField[II + 2][JJ - 1][IR];
    qField[II + 2][JJ + 1][IU] = qField[II + 2][JJ - 1][IU];
    qField[II + 2][JJ + 1][IV] = -qField[II + 2][JJ - 1][IV];
    qField[II + 2][JJ + 1][IP] = qField[II + 2][JJ - 1][IP];

    // 4点， x方向取值和y方向取值的算术平均
    qField[II + 1][JJ + 1][IR] = 0.5 * (qField[II - 1][JJ + 1][IR] + qField[II + 1][JJ - 1][IR]);
    qField[II + 1][JJ + 1][IU] = 0.5 * (qField[II - 1][JJ + 1][IU] + qField[II + 1][JJ - 1][IU]);
    qField[II + 1][JJ + 1][IV] = 0.5 * (qField[II - 1][JJ + 1][IV] + qField[II + 1][JJ - 1][IV]);
    qField[II + 1][JJ + 1][IP] = 0.5 * (qField[II - 1][JJ + 1][IP] + qField[II + 1][JJ - 1][IP]);
}