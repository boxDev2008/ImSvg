#pragma once

#include <imgui.h>
#include <vector>
#include <string>

namespace SVGImGui
{

struct SVGColor { float r=0,g=0,b=0,a=1; bool none=false; };
struct SVGStyle
{
    SVGColor fill       = {0,0,0,1,false};
    SVGColor stroke;
    float strokeWidth   = 1.0f;
    float fillOpacity   = 1.0f;
    float strokeOpacity = 1.0f;
    float opacity       = 1.0f;
    SVGStyle(void) { stroke.none = true; }
};

struct Contour
{
    std::vector<ImVec2> pts;
    bool closed = false;
};
struct Triangle { int a,b,c; };

struct RenderShape
{
    SVGStyle style;
    std::vector<Contour> contours;
    std::vector<std::vector<Triangle>> fillTris;
    bool isLine = false;
    ImVec2 lineA, lineB;
};

struct SVGDocument
{
    float width=100, height=100;
    float viewBox[4]={0,0,0,0};
    bool  hasViewBox=false;
    std::vector<RenderShape> shapes;
};

SVGDocument ParseSVG(const std::string& xml);
SVGDocument LoadSVG(const char* path);

void RenderSVG(ImDrawList* dl, const SVGDocument& doc, ImVec2 origin, float scale = 1.0f, float globalAlpha = 1.0f);
void RenderSVGString(ImDrawList* dl, const std::string& svgXml, ImVec2 origin, float scale = 1.0f, float globalAlpha = 1.0f);

float SvgFitScale(const SVGDocument& doc, float boxW, float boxH);

}