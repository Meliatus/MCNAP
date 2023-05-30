#include "gridTools.h"
#include <glm/glm.hpp>
#include <PGUPV.h>

using namespace PGUPV;

void vdc::sampleVectorField(vdc::Grid<glm::vec2, glm::vec2> &g, std::function<glm::vec2(const glm::vec2 &)> f) {
  for (size_t i = 0; i < g.numSamples(); i++) {
    g.setSampleValue(i, f(g.getSamplePosition(i)));
  }
}

/**

Completa esta función:

La función devuelve un modelo PGUPV con un segmento de línea (dos vértices) por cada muestra de la malla g. 
Dicha línea apunta en la dirección indicada por el valor de la muestra, y tiene una longitud de k.

*/

std::shared_ptr<Model> vdc::computeHedgeHog(const vdc::Grid<glm::vec2, glm::vec2> &g, float k) {
  auto result = std::make_shared<Model>();
  std::vector<glm::vec2> vertices;
  std::vector<glm::vec4> colores;
  std::vector<unsigned int> indices;
  for (size_t i = 0; i < g.numSamples(); i++) {
      glm::vec2 dir = g.getSampleValue(i);
      float mag = glm::length(dir);
      if (mag > 0.0f) {
          dir /= mag;
      }
      auto pos= g.getSamplePosition(i);
      glm::vec2 final = pos + glm::normalize(dir) * k;
      vertices.push_back(pos);
      vertices.push_back(final);
      colores.push_back(glm::vec4(mag, 0.0f, 0.0f, 1.0f));
      colores.push_back(glm::vec4(mag, 0.0f, 0.0f, 1.0f));
      indices.push_back(static_cast<unsigned int>(i * 2));
      indices.push_back(static_cast<unsigned int>(i * 2 + 1));
  }

  auto base = std::make_shared<Mesh>();
  base->addVertices(vertices);
  base->addColors(colores);
  base->addIndices(indices);
  base->addDrawCommand(new DrawElements(GL_LINES, static_cast<GLsizei>(indices.size()), GL_UNSIGNED_INT, 0));
  result->addMesh(base);

  return result;
    
}

/**

Completa esta función:

La función devuelve una malla uniforme del mismo tamaño que la malla de entrada. En vez de ser 
una malla vectorial, será una malla escalar, cuyas muestras serán el valor de la divergencia de 
la muestra correspondiente de la malla de entrada:

*/


std::shared_ptr<vdc::UniformGrid<glm::vec2, float>> vdc::computeDivergence(const vdc::UniformGrid<glm::vec2, glm::vec2> &g) {
    std::vector<int> dims{ g.getNumSamplesPerDimension(0), g.getNumSamplesPerDimension(1) };
    auto result = std::make_shared<vdc::UniformGrid<glm::vec2, float>>(g.getMinCoord(), g.getMaxCoord(), dims);
    for (int i = 0; i < g.numSamples(); i++) {
        int nexty,nextx;
        if ((i + dims[1]) < g.numSamples()) {
            nexty = i + dims[1];
        } else { nexty = i; }
        float distanciay = g.getSamplePosition(nexty).y - g.getSamplePosition(i).y;
        float dy = (g.getSampleValue(nexty).y - g.getSampleValue(i).y) / distanciay;
        if ((i + 1) % dims[1] != 0) {
            nextx = i + 1;
        } else { nextx = i; }
        float distanciax = g.getSamplePosition(nextx).x - g.getSamplePosition(i).x;
        float dx = (g.getSampleValue(nextx).x - g.getSampleValue(i).x) / distanciax;

        float divergencia = dx + dy;
        result->setSampleValue(i, divergencia);
    }

    return result;
}


/*
Completa esta función:

La función devuelve una malla escalar de la misma dimensión que la malla de entrada, donde cada muestra 
contiene la magnitud de la vorticidad de la muestra correspondiente de la malla de entrada.

*/

std::shared_ptr<vdc::UniformGrid<glm::vec2, float>> vdc::computeVorticity(const vdc::UniformGrid<glm::vec2, glm::vec2> &g) {
    std::vector<int> dims{ g.getNumSamplesPerDimension(0), g.getNumSamplesPerDimension(1) };
    auto result = std::make_shared<vdc::UniformGrid<glm::vec2, float>>(g.getMinCoord(), g.getMaxCoord(), dims);
    float dx = (g.getMaxCoord().x - g.getMinCoord().x) / (g.getNumSamplesPerDimension(0) - 1);
    float dy = (g.getMaxCoord().y - g.getMinCoord().y) / (g.getNumSamplesPerDimension(1) - 1);

    for (int i = 0; i < g.numSamples(); i++) {
        int nexty, nextx;
        if ((i + dims[1]) < g.numSamples()) {
            nexty = i + dims[1];
        }
        else { nexty = i; }

        float distanciay = g.getSamplePosition(nexty).y - g.getSamplePosition(i).y;
        float dvx_dy = (g.getSampleValue(nexty).x - g.getSampleValue(i).x) / distanciay;

        if ((i + 1) % dims[1] != 0) {
            nextx = i + 1;
        }
        else { nextx = i; }

        float distanciax = g.getSamplePosition(nextx).x - g.getSamplePosition(i).x;
        float dvy_dx = (g.getSampleValue(nextx).y - g.getSampleValue(i).y) / distanciax;

        float vorticidad = dvy_dx - dvx_dy;
        result->setSampleValue(i, vorticidad);
    }
    return result;
}


/*

Completa la siguiente función. 

La función devuelve una malla de PGUPV con los vértices, colores y draw command 
necesario para dibujar la línea de corriente que empieza en p0. El paso de integración 
se pasa en el parámetro dt, maxT es el tiempo máximo de integración y maxL es 
la longitud máxima de la línea de corriente.Por último, la línea se dibujará 
del color indicado por el último parámetro.
*/

std::shared_ptr<PGUPV::Mesh> vdc::computeStreamline(const vdc::UniformGrid<glm::vec2, glm::vec2> &g, glm::vec2 &p0, float dt, float maxT, float maxL, glm::vec4 color) {
    auto result = std::make_shared<PGUPV::Mesh>();

    /*glm::vec2 v, p, q;
    int c = g.findCell(p0); // celda inicial
    float t = 0, l = 0;
    int i = 0;
    while (true) {
        if (c == -1 || t > maxT || l > maxL) break; // criterios de parada
        result->addVertices({ p0.x, p0.y }); // añadir vértice a la línea
        if (i > 0) {
            result->addDrawCommand(PGUPV::DrawCommand::LINE_STRIP, i - 1, 2); // añadir nuevo comando de dibujo
            result->addColor(color, i - 1); // añadir color a la línea de corriente
            l += glm::length(p - p0); // actualizar longitud
        }
        p = p0;
        g.interpolateC1Square(p0, v);
        p0 += v * dt; // Integración de Euler
        g.world2cell(p0, q);
        if (q.x < 0 || q.x > 1 || q.y < 0 || q.y > 1) // ¿nos hemos salido de la celda?
            c = g.findCell(p0); // encontrar la nueva celda
        t += dt;
        i++;
    }*/

    return result;
}