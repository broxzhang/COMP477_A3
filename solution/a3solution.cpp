#include "a3solution.h"

#include "dependencies/Eigen/Dense"
#include "QDebug"
#include "QElapsedTimer"

using Eigen::Vector2f;
using Eigen::Vector4f;
using Eigen::VectorXf;
using Eigen::MatrixXd;
using Eigen::MatrixXf;
using std::vector;

A3Solution::A3Solution(std::vector<Joint2D*>& joints, std::vector<Spring2D*>& springs, float& gravity, float& positional_damping, float& mass, float& timestep, bool& implicit, float& stiffness)
    :m_joints(joints),
    m_links(springs),
    m_gravity(gravity),
    m_positional_damping(positional_damping),
    m_mass(mass),
    m_timestep(timestep),
    m_implicit(implicit),
    m_stiffness(stiffness),
    m_JointsStates()
    {
    m_additionalSprings vector<Spring2D*>();
}


void A3Solution::update(Joint2D* selected, QVector2D mouse_pos){
    selected->set_position(mouse_pos);

    m_JointsStates[selected].x = Vector2f(mouse_pos.x(), mouse_pos.y());
    m_JointsStates[selected].v = Vector2f(0.0f, 0.0f);
}

void A3Solution::update(){
    // Students implement solution here

    for (Joint2D* j : m_joints) {
        if (m_JointsStates.find(j) == m_JointsStates.end()) {
            y initialState;
            QVector2D jointPosition = j->get_position();
            initialState.x = Vector2f(jointPosition.x(), jointPosition.y());
            initialState.v = Vector2f(0,0f, 0.0f);

            m_JointsStates[j] = initialState;
        }
    }

    for (auto it = m_JointsState.begin(); it != m_JointsState.end(); it++) {
        if (std::find(m_joints.begin(), m_joints.end(), it->first) == m_joints.end()) {
            m_JointsState.erase(it);
        } else {
            Joint2D* j = it->first;
            if(m_implicit){
                m_JointsStates[j] = SemiImplicitEuler(j);
            } else {
                m_JointsStates[j] = ExplicitEuler(j);
            }
        }
    }

    for (Joint2D* j : m_joints) {
        j->set_position(m_JointsStates[j].x);
    }
}

A3Solution::y A3Solution::ExplicitEuler(Joint2D* j) {
    if (!j->is_locked()) {
        y k = m_JointsStates[j];
        y k1;
        k1.x = k.x + m_timestep * k.v;
        k1.v = k.v + m_timestep * (total_force(j) / m_mass);

        return k1;
    }
    return m_JointsStates[j];
}

A3Solution::y A3Solution::SemiImplicitEuler(Joint2D* j) {
    if (!j->is_locked()) {
        y k = m_JointsStates[j];
        y k1;
        k1.v = k.v + m_timestep * (total_force(j) / m_mass);
        k1.x = k.x + m_timestep * k1.v;

        return k1;
    }
    return m_JointsStates[j];
}

Vector2f A3Solution::total_force(Joint2D* j) {
    Vector2f g(0.0f, m_gravity * m_mass);
    Vector2f spring(0.0f, 0.0f);
    Vector2f damping(-m_positional_damping * m_JointsStates[j].v);

    for (Spring2D* s : j->get_springs()) {

        QVector2D j1 = j->get_position();
        QVector2D j2 = s->get_other_joint(j)->get_position();

        Vector2f x1(j1.x(), j1.y());
        Vector2f x2(j2.x(), j2.y());

        Vector2f x = x1 - x2;
        float delta_norm = x.norm();

        spring +=  -m_stiffness * (x / (delta_norm * delta_norm) * (delta_norm - s->get_rest_length()));
    }
}

void A3Solution::test_eigen_library(){

    // create a simple matrix 5 by 6
    MatrixXd mat(5,6);

    // Fills in matrix
    // Important Note: Eigen matrices are row major
    // so mat(0,1) references the 0-th column and 1-th row
    for(unsigned int row=0;row<mat.rows();row++){
        for(unsigned int col=0;col<mat.cols();col++){
            mat(row,col) = row+col;
        }
    }

    // create the pseudoinverse
    MatrixXd pseudo_inv = mat.completeOrthogonalDecomposition().pseudoInverse();

    // print the pseudoinverse
    std::cout << "--------------------------" << std::endl;
    std::cout << pseudo_inv << std::endl;
    std::cout << "--------------------------" << std::endl;

}
